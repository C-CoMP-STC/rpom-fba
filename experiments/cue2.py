import os
import pickle
import warnings
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
from cobra.io import load_model, read_sbml_model
from matplotlib import pyplot as plt

from experiments.cue import plot_bge, plot_result
from experiments.experiment import Experiment
from experiments.fast_dFBA import (ConstantBounds, dFBA, make_bge_listener,
                                   make_boundary_listener, make_cue_listener,
                                   make_growth_rate_listener,
                                   make_shadow_price_listener,
                                   plot_shadow_prices, setup_drawdown)
from parameters.drawdown import *
from utils.cobra_utils import get_or_create_exchange
from utils.units import u

C_PER_GLUCOSE = 6
C_PER_ACETATE = 2

OUTDIR = "out/dFBA/CUE2"


class CUE_Experiment_2(Experiment):
    def __init__(
        self, model, biomass_id, ex_glc="EX_glc", ex_ace="EX_ace", out=OUTDIR
    ) -> None:
        self.model = model
        self.biomass_id = biomass_id
        self.ex_glc = model.reactions.get_by_id(ex_glc)
        self.ex_ace = model.reactions.get_by_id(ex_ace)
        self.dt = 0.1

        # Ensure output directory exists
        self.out = out
        os.makedirs(out, exist_ok=True)
        self.data_out = os.path.join(
            out, datetime.now().strftime("%d-%b-%Y_%H:%M:%S"))
        os.makedirs(self.data_out, exist_ok=True)

    def unpack_condition_data(self, condition, condition_data):
        # Get initial conditions
        initial_glucose, initial_acetate = tuple(
            map(lambda x: x.magnitude, condition))

        # Get biomass timeseries from data
        b_t, b = condition_data["mean"]["b_t"], condition_data["mean"]["b_s"]
        initial_biomass = b[0]
        tmax = b_t.max()

        initial = np.array([initial_biomass, initial_glucose, initial_acetate])

        # Get glucose and acetate traces from data
        substrate_trace = {}
        for s0, s_t, s in [
            (initial_glucose, "g_t", "g_s"),
            (initial_acetate, "a_t", "a_s")
        ]:
            substrate_trace[s_t] = (
                condition_data["mean"][s_t] if s0 != 0 else np.array([0, tmax])
            )
            substrate_trace[s] = (
                condition_data["mean"][s] if s0 != 0 else np.array([0, 0])
            )

        return (
            initial,
            tmax,
            b_t,
            b,
            substrate_trace["g_t"],
            substrate_trace["g_s"],
            substrate_trace["a_t"],
            substrate_trace["a_s"],
        )

    def run_condition(self, condition, condition_data, save_data=True):
        initial, tmax, b_t, b, g_t, g_s, a_t, a_s = self.unpack_condition_data(
            condition, condition_data
        )
        initial_biomass, initial_glucose, initial_acetate = initial

        dynamic_medium = [
            ConstantBounds(self.ex_glc, "Glucose[e]", -4),
            ConstantBounds(self.ex_ace, "ACET[e]", -2)
        ]

        # Override initial biomass and dynamic medium if given
        initial_biomass = condition_data.get(
            "initial_biomass", initial_biomass)
        dynamic_medium = condition_data.get("dynamic_medium", dynamic_medium)
        for bounds in dynamic_medium:
            bounds.dt = self.dt

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t, y, l = dFBA(
                self.model,
                self.biomass_id,
                initial_biomass,
                np.array([initial_glucose, initial_acetate]),
                dynamic_medium,
                tmax,
                terminate_on_infeasible=False,
                listeners=[
                    #    make_shadow_price_listener(
                    #        self.model, ["Glucose[e]", "ACET[e]"], dynamic_medium),
                    # TODO: bge listener causes freezing on some conditions
                    # make_bge_listener(
                    # self.model, self.biomass_id, dynamic_medium),
                    # make_growth_rate_listener(
                    #     self.model, self.biomass_id, dynamic_medium
                    # ),
                    # make_boundary_listener(
                    #     self.model, self.biomass_id, dynamic_medium
                    # )
                ],
                dt=self.dt,
                integrator="runge_kutta",
                use_cache=True
            )

        if save_data:
            self.save_data(
                os.path.join(
                    self.data_out,
                    f"data_{initial_glucose:.2f}mM_glucose_{initial_acetate:.2f}mM_acetate.csv",
                ),
                t, y, l,
            )

        return t, y, l

    def save_data(self, filepath, t, y, l):
        save_data = pd.DataFrame(
            {
                "time (h)": t,
                "Biomass (g/L)": y[:, 0],
                "Glucose (mM)": y[:, 1],
                "Acetate (mM)": y[:, 2],
            }
        )
        if len(l) > 0:
            pass
        #     sp = pd.concat(l[0], axis=1, ignore_index=True).T
        #     save_data = pd.concat((save_data, sp), axis=1)
        save_data.to_csv(filepath, index=False)


def main():
    BIOMASS_ID = "Rpom_hwa_biomass"
    DATA_FILE = "data/clean/CUE2/dFBA.pkl"

    # Load and set up model
    model = read_sbml_model("model/Rpom_05.xml")
    setup_drawdown(model)
    ex_ace = get_or_create_exchange(model, "ACET[e]")

    # Turn on maintenance
    atpm = model.reactions.get_by_id("ATPM")
    atpm.bounds = (30, 30)

    # Replace all ubiquinones with ubiquinone-8, same with ubiquinols. TODO: curate this more, use whichever R. pom mostly uses.
    ubiquinones = [model.metabolites.get_by_id(ubi)
                   for ubi in [
        "UBIQUINONE-6[c]",
        "UBIQUINONE-8[c]",
        "UBIQUINONE-9[c]",
        "UBIQUINONE-10[c]"]]
    ubiquinols = [model.metabolites.get_by_id(ubi)
                  for ubi in [
        "UBIQUINOL-30[c]",
        "CPD-9956[c]",
        "CPD-9957[c]",
        "CPD-9958[c]"
    ]]
    ubiquinone8 = model.metabolites.get_by_id("UBIQUINONE-8[c]")
    ubiquinol8 = model.metabolites.get_by_id("CPD-9956[c]")
    for rxn in model.reactions:
        metabolites_copy = rxn.metabolites.copy()
        for met, coeff in metabolites_copy.items():
            if met in ubiquinones:
                rxn.add_metabolites({met: -coeff, ubiquinone8: coeff})
            if met in ubiquinols:
                rxn.add_metabolites({met: -coeff, ubiquinol8: coeff})

    cue_experiment = CUE_Experiment_2(model, BIOMASS_ID, ex_ace=ex_ace.id)
    # cue_experiment = CUE_Experiment(model, "BIOMASS_Ec_iJO1366_core_53p95M", ex_ace="EX_ac_e", ex_glc="EX_glc__D_e")

    # Load condition data
    with open(DATA_FILE, "rb") as f:
        data = pickle.load(f)
    cue_experiment.conditions = data

    # Use fitted parameters and initial conditions
    # {(<Quantity(2.0, 'millimolar')>, <Quantity(0.0, 'millimolar')>): ((3.857788639318984e-06, -3.7498565796153476, -20.0), 4.864342853437677)}
    # {(<Quantity(0.0, 'millimolar')>, <Quantity(6.0, 'millimolar')>): ((1.7746685352765677e-06, -4.0, -20.65117661222312), 4.604088166504376)}

    fitted_params = {
        (2.0, 0.0): np.array([3.857788639318984e-06, -5, -0]),
        (0.0, 6.0): np.array([5e-05, -0, -18]),
        (2/3, 4.0): np.array([1.852104325893396e-05, -2.5 * 2/3, -3 * 4.]),
        (4/3, 2.0): np.array([1.5705806351090746e-05, -2.5 * 4/3, -3 * 2.]),
    }
    for k, v in fitted_params.items():
        key = [(g, a) for (g, a) in cue_experiment.conditions if g.magnitude ==
               k[0] and a.magnitude == k[1]][0]
        cue_experiment.conditions[key]["initial_biomass"] = v[0]
        cue_experiment.conditions[key]["dynamic_medium"] = [
            ConstantBounds(cue_experiment.ex_glc, "Glucose[e]", v[1]),
            ConstantBounds(cue_experiment.ex_ace, "ACET[e]", v[2])
        ]

    # Run dFBA
    # cue_experiment.dt = 1
    # TODO: delete
    # cue_experiment.conditions = {
    #     (g, a): v
    #     for (g, a), v in cue_experiment.conditions.items()
    #     if g == 0 and a.magnitude > 0
    # }
    results = cue_experiment.run()

    for condition, (t, y, l) in results.items():
        # Plot shadow prices
        # fig, axs = plot_shadow_prices(l[0], t)
        # axs[0].set_title(
        #     f"{initial_glucose.magnitude:.2f}mM Glucose, {initial_acetate.magnitude:.2f} mM Acetate")
        # fig.set_size_inches(5, len(axs))
        # fig.savefig(os.path.join(
        #     OUTDIR, f"{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_shadow_prices.png"))

        # Plot data
        initial_glucose, initial_acetate = condition
        condition_data = cue_experiment.conditions[condition]
        fig, (ax, _) = plot_result(t, y, condition_data, mass_units=False)
        ax.set_yscale("log")

        # plot bge
        # ax3 = ax.twinx()
        # bge, mu, boundary = l
        # plot_bge(t, bge, ax=ax3)

        # Save figure
        fig.set_size_inches(6, 3)
        fig.tight_layout()
        fig.savefig(
            os.path.join(
                OUTDIR,
                f"{model.id}_{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_dFBA.png",
            )
        )

        # Plot mu
        # fig, ax = plt.subplots()
        # b = y[:, 0]
        # ax.plot(t, b, "b")
        # integ = [b[0]]
        # dt = t[1] - t[0]
        # for mu in mu:
        #     integ.append(integ[-1] + mu * integ[-1] * dt)
        # ax.plot(t, integ[:-1], "b--")

        # ax.set_xlabel("Biomass (g/L)")
        # fig.tight_layout()
        # fig.savefig(os.path.join(
        #     OUTDIR,
        #     f"{model.id}_mu_test_{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_dFBA.png",
        # ))

        # # Make boundary figure
        # traces = defaultdict(list)
        # for timepoint in boundary:
        #     for reaction, flux in timepoint:
        #         traces[reaction.id].append(flux)

        # n_traces = len(traces)

        # if n_traces > 0:
        #     fig, axs = plt.subplots(n_traces, 1)
        #     for ax, (reaction, trace) in zip(axs, traces.items()):
        #         ax.plot(trace)
        #         ax.set_ylabel(reaction)

        #     fig.set_size_inches(6, n_traces * 2)
        #     fig.tight_layout()
        #     fig.savefig(os.path.join(OUTDIR, f"{model.id}_boundary_{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_dFBA.png"))


if __name__ == "__main__":
    main()
