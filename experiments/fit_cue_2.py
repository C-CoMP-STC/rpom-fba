import pickle

import numpy as np
from cobra.io import read_sbml_model
from experiments.cue2 import CUE_Experiment_2
from experiments.fast_dFBA import ConstantBounds, setup_drawdown
from utils.cobra_utils import get_or_create_exchange
from utils.units import u
from utils.math import get_interpolator


def main():
    BIOMASS_ID = "Rpom_hwa_biomass"
    DATA_FILE = "data/clean/CUE2/dFBA.pkl"

    # Load and set up model
    model = read_sbml_model("model/Rpom_05.xml")
    setup_drawdown(model)
    ex_ace = get_or_create_exchange(model, "ACET[e]")

    cue_experiment = CUE_Experiment_2(model, BIOMASS_ID, ex_ace=ex_ace.id)

    # Load condition data
    with open(DATA_FILE, "rb") as f:
        condition_data = pickle.load(f)

    def initialize_model(log_b0, dgdt, dadt):
        # Assuming only one condition at a time
        k = list(cue_experiment.conditions.keys())[0]
        cue_experiment.conditions[k]["initial_biomass"] = np.exp(log_b0)
        cue_experiment.conditions[k]["dynamic_medium"] = [
            ConstantBounds(cue_experiment.ex_glc, "Glucose[e]", dgdt),
            ConstantBounds(cue_experiment.ex_ace, "ACET[e]", dadt),
        ]
        return cue_experiment

    result = {}
    for k, data in condition_data.items():
        glucose, acetate = k
        print(f"Condition: {k} =====================")

        if not (glucose > 0 and acetate > 0):
            continue

        # only run one condition at a time
        cue_experiment.conditions = {k: data}

        def sse(t, b, g=None, a=None):
            b_t_true = data["mean"]["b_t"]
            b_true = data["mean"]["b_s"]
            g_t_true = data["mean"]["g_t"]
            g_true = np.nan_to_num(data["mean"]["g_s"])
            a_t_true = data["mean"]["a_t"]
            a_true = np.nan_to_num(data["mean"]["a_s"])

            b_t = get_interpolator(t, b)
            b_hat = b_t(b_t_true)

            if g is not None:
                g_t = get_interpolator(t, g)
                g_hat = g_t(g_t_true)
            else:
                g_hat = np.zeros_like(g_true)

            if a is not None:
                a_t = get_interpolator(t, a)
                a_hat = a_t(a_t_true)
            else:
                a_hat = np.zeros_like(a_true)

            # weight by 10?? so that b and g, a errors are on the same scale
            return (
                ((np.log(b_hat) - np.log(b_true)) ** 2).sum()
                + ((g_hat - g_true) ** 2).sum()
                + ((a_hat - a_true) ** 2).sum()
            )

        def grad(
            log_b0,
            dgdt,
            dadt,
            steps=(0.01, 0.01, 0.01),
            current_error=None,
            skip_g=False,
            skip_a=False,
        ):
            if current_error is None:
                (t, y, _) = initialize_model(log_b0, dgdt, dadt).run()[k]
                b = y.T[0]
                g = y.T[1]
                a = y.T[2]
                f = sse(t, b, g, a)
            else:
                f = current_error

            grad = np.array([0, 0, 0])
            for i, step in enumerate(steps):
                if skip_g and i == 1:
                    continue
                if skip_a and i == 2:
                    continue
                params = [log_b0, dgdt, dadt]
                params[i] += step

                (t, y, _) = initialize_model(*params).run()[k]
                b = y.T[0]
                g = y.T[1]
                a = y.T[2]
                grad[i] = (sse(t, b, g, a) - f) / step

            print(grad)
            return grad

        # Include parameters in condition dictionary
        # parameterize with log_b0 to match scale of data
        log_b0 = np.log(data["mean"]["b_s"][0])
        dgdt = -4
        dadt = -2

        # Get error of initial guess
        initialize_model(log_b0, dgdt, dadt)
        cue_experiment.dt = 1
        (t, y, _) = cue_experiment.run()[k]
        b = y.T[0]
        g = y.T[1]
        a = y.T[2]
        error = sse(t, b, g, a)
        print(f"Initial: {(np.exp(log_b0), dgdt, dadt)} ({error=})")

        model = (log_b0, dgdt, dadt)
        best_model = model
        best_error = error

        iteration = 0
        MAXITER = 50
        dt_schedule = [1, 0.5, 0.1]
        convergence_target = [1, 0.1, 0.01]
        step_schedule = np.logspace(-3, -4, num=len(dt_schedule) * MAXITER)  # -5?

        # Use backtracking line search instead to adjust step size
        # alpha = 0.5
        # beta = 0.5
        # step = step_schedule[iteration]
        # g = grad(
        #     *model,
        #     current_error=error,
        #     skip_a=(acetate == 0),
        #     skip_g=(glucose == 0),
        # )
        # while True:
        #     new_model = model - step * g
        #     (t, y, _) = initialize_model(*new_model).run()[k]
        #     b = y.T[0]
        #     g = y.T[1]
        #     a = y.T[2]
        #     new_error = sse(t, b, g, a)
        #     if new_error < error + alpha * step * np.dot(g, g):
        #         break
        #     step *= beta
        #     model = new_model
        #     error = new_error


        for dt, target in zip(dt_schedule, convergence_target):
            print(f"Trying dt = {dt}")
            cue_experiment.dt = dt

            diff = 1e10
            while abs(diff) > target and (iteration / MAXITER) < 1:
                g = grad(
                    *model,
                    current_error=error,
                    skip_a=(acetate == 0),
                    skip_g=(glucose == 0),
                )

                model -= step_schedule[iteration] * g

                (t, y, _) = initialize_model(*model).run()[k]
                b = y.T[0]
                g = y.T[1]
                a = y.T[2]
                new_error = sse(t, b, g, a)

                diff = new_error - error
                error = new_error
                print(
                    f"({iteration}) Model: {(np.exp(model[0]), model[1], model[2])}, {error=} ({diff=})"
                )

                if error < best_error:
                    best_error = error
                    best_model = model

                iteration += 1
        
        result[k] = ((np.exp(best_model[0]), best_model[1], best_model[2]), best_error)
        print(
            f"Best model for {k}: {(np.exp(best_model[0]), best_model[1], best_model[2])} ({best_error})"
        )

    print("FINAL RESULTS =====================")
    print(result)


if __name__ == "__main__":
    main()
