from dataclasses import dataclass
import pickle
import multiprocessing
import os

import numpy as np
from cobra.io import read_sbml_model
from experiments.cue2 import CUE_Experiment_2
from experiments.fast_dFBA import ConstantBounds, setup_drawdown
from utils.cobra_utils import get_or_create_exchange
from utils.units import u
from utils.math import get_interpolator

global condition_data


def initialize_model(cue_experiment, log_b0, dgdt, dadt):
    # Assuming only one condition at a time
    k = list(cue_experiment.conditions.keys())[0]
    cue_experiment.conditions[k]["initial_biomass"] = np.exp(log_b0)
    cue_experiment.conditions[k]["dynamic_medium"] = [
        ConstantBounds(cue_experiment.ex_glc, "Glucose[e]", dgdt),
        ConstantBounds(cue_experiment.ex_ace, "ACET[e]", dadt),
    ]
    return cue_experiment


def sse(t, b, g, a, condition):
    global condition_data

    data = condition_data[condition]

    b_t_true = data["mean"]["b_t"]
    b_true = data["mean"]["b_s"]
    b_hat = get_interpolator(t, b)(b_t_true)

    g_t_true = data["mean"]["g_t"]
    g_true = np.nan_to_num(data["mean"]["g_s"])
    g_hat = get_interpolator(t, g)(g_t_true)

    a_t_true = data["mean"]["a_t"]
    a_true = np.nan_to_num(data["mean"]["a_s"])
    a_hat = get_interpolator(t, a)(a_t_true)

    return (
        np.nan_to_num((b_hat - b_true) ** 2).sum()
        + np.nan_to_num((g_hat - g_true) ** 2).sum()
        + np.nan_to_num((a_hat - a_true) ** 2).sum()
    )


@dataclass
class GradStepParams:
    i: int
    step: float
    cue_experiment: CUE_Experiment_2
    skip_g: bool
    skip_a: bool
    log_b0: float
    dgdt: float
    dadt: float
    condition: tuple
    data: dict
    current_error: float


def do_grad_step(grad_params: GradStepParams):
    """Helper function to calculate gradient in parallel"""
    if grad_params.skip_g and grad_params.i == 1:
        return 0
    if grad_params.skip_a and grad_params.i == 2:
        return 0
    params = [grad_params.log_b0, grad_params.dgdt, grad_params.dadt]
    params[grad_params.i] += grad_params.step

    (t, y, _) = initialize_model(
        grad_params.cue_experiment,
        *params
    ).run_condition(grad_params.condition, grad_params.data, save_data=False)

    b = y.T[0]
    g = y.T[1]
    a = y.T[2]
    return (sse(t, b, g, a, grad_params.condition) - grad_params.current_error) / grad_params.step


def fit_condition(x, biomass_id="Rpom_hwa_biomass", ex_ace_id="EX_ac", parallel_gradient=False):
    model, condition, data = x
    glucose, acetate = condition

    print(f"[{os.getpid()}] Condition: {condition} =====================")

    # only run one condition at a time
    cue_experiment = CUE_Experiment_2(model, biomass_id, ex_ace=ex_ace_id)
    cue_experiment.conditions = {condition: data}

    def grad(
        log_b0,
        dgdt,
        dadt,
        steps=(1e-2, 1e-2, 1e-2),
        current_error=None,
        skip_g=False,
        skip_a=False,
        parallel=False,
    ):
        if current_error is None:
            (t, y, _) = initialize_model(cue_experiment, log_b0, dgdt,
                                         dadt).run_condition(condition, data, save_data=False)
            b = y.T[0]
            g = y.T[1]
            a = y.T[2]
            f = sse(t, b, g, a, condition)
        else:
            f = current_error

        grad = np.array([0., 0., 0.])
        if parallel:
            # def do_grad_step(x):
            #     i, step = x
            #     if skip_g and i == 1:
            #         return 0
            #     if skip_a and i == 2:
            #         return 0
            #     params = [log_b0, dgdt, dadt]
            #     params[i] += step

            #     (t, y, _) = initialize_model(cue_experiment, *params).run_condition(condition, data, save_data=False)
            #     b = y.T[0]
            #     g = y.T[1]
            #     a = y.T[2]
            #     return (sse(t, b, g, a, condition) - f) / step

            grad_params = [GradStepParams(i, step, cue_experiment, skip_g, skip_a, log_b0, dgdt, dadt, condition, data, f)
                           for i, step in enumerate(steps)]
            with multiprocessing.Pool() as pool:
                grad = np.array(pool.map(do_grad_step, grad_params))
        else:
            for i, step in enumerate(steps):
                if skip_g and i == 1:
                    continue
                if skip_a and i == 2:
                    continue
                params = [log_b0, dgdt, dadt]
                params[i] += step

                (t, y, _) = initialize_model(cue_experiment, *
                                             params).run_condition(condition, data, save_data=False)
                b = y.T[0]
                g = y.T[1]
                a = y.T[2]
                grad[i] = (sse(t, b, g, a, condition) - f) / step

        print(grad)
        return grad

    # Include parameters in condition dictionary
    # parameterize with log_b0 to match scale of data
    log_b0 = np.log(data["mean"]["b_s"][0])
    dgdt = -4
    dadt = -20

    # Get error of initial guess
    initialize_model(cue_experiment, log_b0, dgdt, dadt)
    cue_experiment.dt = 1
    (t, y, _) = cue_experiment.run_condition(condition, data, save_data=False)
    # (t, y, _) = cue_experiment.run()[condition]
    b = y.T[0]
    g = y.T[1]
    a = y.T[2]
    error = sse(t, b, g, a, condition)
    print(f"[{os.getpid()}] Initial: {(np.exp(log_b0), dgdt, dadt)} ({error=})")

    model = (log_b0, dgdt, dadt)
    best_model = model
    best_error = error

    iteration = 0
    MAXITER = 30
    # dt_schedule = [1, 0.5, 0.1]
    dt_schedule = [1, 0.5]
    convergence_target = [1, 0.01]
    step_schedule = np.logspace(-2, -3, num=len(dt_schedule) * MAXITER)

    for i, (dt, target) in enumerate(zip(dt_schedule, convergence_target)):
        print(f"[{os.getpid()}] Trying dt = {dt}")
        cue_experiment.dt = dt

        diff = float("inf")
        while abs(diff) > target and (iteration / MAXITER) < (i + 1):
            g = grad(
                *model,
                current_error=error,
                skip_a=(acetate == 0),
                skip_g=(glucose == 0),
                parallel=parallel_gradient
            )

            model -= step_schedule[iteration] * g

            (t, y, _) = initialize_model(cue_experiment, *
                                         model).run_condition(condition, data, save_data=False)
            b = y.T[0]
            g = y.T[1]
            a = y.T[2]
            new_error = sse(t, b, g, a, condition)

            diff = new_error - error
            error = new_error
            print(
                f"[{os.getpid()}] ({iteration}) Model: {(np.exp(model[0]), model[1], model[2])}, {error=} ({diff=})"
            )

            if error < best_error:
                best_error = error
                best_model = model

            iteration += 1

    best_model = (np.exp(best_model[0]), best_model[1], best_model[2])
    return best_model, best_error


def main():
    global condition_data

    # Flags to toggle parallelization of conditions and gradient calculation
    # Currently, only one of these can be set to True at a time
    PARALLEL_CONDITIONS = False
    PARALLEL_GRADIENT = True

    DATA_FILE = "data/clean/CUE2/dFBA.pkl"

    # Load and set up model
    model = read_sbml_model("model/Rpom_05.xml")
    setup_drawdown(model)
    ex_ace = get_or_create_exchange(model, "ACET[e]")

    # Turn on ATPM
    atpm = model.reactions.get_by_id("ATPM")
    atpm.bounds = (30, 30)

    # Load condition data
    with open(DATA_FILE, "rb") as f:
        condition_data = pickle.load(f)

    # TODO: remove
    condition_data = {
        (g, a):  v
        for (g, a), v in condition_data.items()
        if 0 == g.magnitude and a > 0
    }

    # Replace zeros with lowest non-zero value
    for k, v in condition_data.items():
        b_s = v["raw"]["raw_b"]
        v["raw"]["raw_b"] = np.where(b_s == 0, b_s[b_s > 0].min(), b_s)

    if PARALLEL_CONDITIONS:
        with multiprocessing.Pool() as pool:
            result = pool.map(
                fit_condition,
                [(model, k, v) for k, v in condition_data.items()])
    else:
        result = [fit_condition((model, k, v), parallel_gradient=PARALLEL_GRADIENT)
                  for k, v in condition_data.items()]

    result = dict(zip(condition_data.keys(), result))

    print("FINAL RESULTS =====================")
    print(result)


if __name__ == "__main__":
    main()
