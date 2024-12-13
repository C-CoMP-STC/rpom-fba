from collections import defaultdict
import os
import pickle
from datetime import datetime

import numpy as np
import pandas as pd
from cobra.io import read_sbml_model, load_model
from matplotlib import pyplot as plt

from experiments.experiment import Experiment
from experiments.fast_dFBA import (
    ConstantBounds,
    dFBA,
    make_bge_listener,
    make_cue_listener,
    make_shadow_price_listener,
    plot_shadow_prices,
    make_boundary_listener,
    make_growth_rate_listener
)
from parameters.drawdown import *
from utils.units import u
from utils.cobra_utils import get_or_create_exchange


C_PER_GLUCOSE = 6
C_PER_ACETATE = 2

OUTDIR = "out/dFBA/CUE"


class CUE_Experiment(Experiment):
    def __init__(
        self, model, biomass_id, ex_glc="EX_glc", ex_ace="EX_ace", out=OUTDIR
    ) -> None:
        self.model = model
        self.biomass_id = biomass_id
        self.ex_glc = model.reactions.get_by_id(ex_glc)
        self.ex_ace = model.reactions.get_by_id(ex_ace)

        # Ensure output directory exists
        self.out = out
        os.makedirs(out, exist_ok=True)
        self.data_out = os.path.join(
            out, datetime.now().strftime("%d-%b-%Y_%H:%M:%S"))
        os.makedirs(self.data_out)

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
            ConstantBounds(self.ex_glc, "Glucose[e]", -.75),
            ConstantBounds(self.ex_ace, "ACET[e]", -.75),
            # BoundFromData(self.ex_glc, "Glucose[e]", g_t, g_s, b_t, b),
            # BoundFromData(self.ex_ace, "ACET[e]", a_t, a_s, b_t, b),
        ]

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
                # make_bge_listener(
                    # self.model, self.biomass_id, dynamic_medium),
                make_bge_listener(
                    self.model, self.biomass_id, dynamic_medium
                )
                # make_growth_rate_listener(
                #     self.model, self.biomass_id, dynamic_medium
                # ),
                # make_boundary_listener(
                #     self.model, self.biomass_id, dynamic_medium
                # )
            ],
            dt=0.1,
            integrator="runge_kutta",
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


def plot_mean_min_max(
    t, data, across_axis=0, color=None, linestyle="--", alpha=0.2, ax=None
):
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(t, data.mean(axis=across_axis), color=color, linestyle=linestyle)
    ax.fill_between(
        t,
        data.min(axis=across_axis),
        data.max(axis=across_axis),
        color=color,
        alpha=alpha,
    )
    return ax


def plot_result(t, y, condition_data, mass_units=True):
    glc_mw = 180.15588 * u.g / u.mol
    ace_mw = 59.04402 * u.g / u.mol

    fig, ax = plt.subplots()

    # Plot simulation results
    ax.plot(t, y[:, 0], color="b", label="Biomass")
    ax2 = ax if mass_units else plt.twinx(ax)

    if mass_units:
        ax2.plot(
            t,
            ((y[:, 1] * u.mM) * glc_mw).to("g/L").magnitude,
            color="r",
            label=f"Glucose (simulation)",
        )
        ax2.plot(
            t,
            ((y[:, 2] * u.mM) * ace_mw).to("g/L").magnitude,
            color="orange",
            label=f"Acetate (simulation)",
        )
    else:
        ax2.plot(t, y[:, 1], color="r", label=f"Glucose (simulation)")
        ax2.plot(t, y[:, 2], color="orange", label=f"Acetate (simulation)")

    # Plot biomass from data
    mass_t = condition_data["mean"]["b_t"]
    mass_raw = condition_data["raw"]["raw_b"]
    plot_mean_min_max(mass_t, mass_raw, color="b", ax=ax)

    # Plot substrates from data
    glucose_t = condition_data["raw"]["raw_g_t"]
    glucose = condition_data["raw"]["raw_g_s"] * (
        1 if not mass_units else (u.mM * glc_mw).to("g/L").magnitude
    )
    acetate_t = condition_data["raw"]["raw_a_t"]
    acetate = condition_data["raw"]["raw_a_s"] * (
        1 if not mass_units else (u.mM * ace_mw).to("g/L").magnitude
    )

    for t, substrate, color in [
        (glucose_t, glucose, "r"),
        (acetate_t, acetate, "orange"),
    ]:
        if substrate.size > 0:
            plot_mean_min_max(t, substrate, color=color, ax=ax2)

    # Plot BGE
    ax3 = ax.twinx()
    bge_t = condition_data["raw"]["raw_bge_t"]
    bge = condition_data["raw"]["raw_bge"]
    plot_mean_min_max(bge_t, bge, color="k", ax=ax3)

    ax.set_ylabel("Biomass (g/L)", color="b")
    ax2.set_ylabel(f"Substrate ({'g/L' if mass_units else 'mM'})", color="r")
    ax2.legend()

    return fig, [ax, ax2]


def plot_bge(t, c, ax=None):
    if ax is None:
        _, ax = plt.subplots()

    # Plot bge
    ax.plot(t, c, "k--")
    ax.set_ylabel("BGE")
    ax.set_ylim(0, 1)
    ax.spines["right"].set_position(("outward", 50))

    return ax


def main():
    BIOMASS_ID = "Rpom_hwa_biomass"
    DATA_FILE = "data/clean/CUE/dFBA.pkl"

    # Load and set up model
    model = read_sbml_model("model/Rpom_05.xml")
    ex_ace = get_or_create_exchange(model, "ACET[e]")
    # model = load_model("iJO1366")

    cue_experiment = CUE_Experiment(model, BIOMASS_ID, ex_ace=ex_ace.id)
    # cue_experiment = CUE_Experiment(model, "BIOMASS_Ec_iJO1366_core_53p95M", ex_ace="EX_ac_e", ex_glc="EX_glc__D_e")

    # Load condition data
    with open(DATA_FILE, "rb") as f:
        data = pickle.load(f)
    cue_experiment.conditions = data

    # Run dFBA
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

        # plot CUE (temp)
        # TODO: delete
        ax3 = ax.twinx()
        cue = l[0]
        plot_bge(t, cue, ax=ax3)


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
