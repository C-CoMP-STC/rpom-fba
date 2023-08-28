import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import matplotlib

from cobra.io import read_sbml_model
from experiments.fast_dFBA import dFBA, setup_drawdown, MichaelisMentenBounds

from utils.cobra_utils import get_or_create_exchange
from parameters.drawdown import *

matplotlib.use("Agg")


EXPERIMENT_VOLUME = 0.1  # L
# TODO: Estimated from E coli, get a better number if possible and move to a single source
MASS_PER_CELL = 0.95  # pg
C_PER_GLUCOSE = 6
C_PER_ACETATE = 2


def plot_result(t, y, initial, V_max, t_max, data):
    condition_data = data[(data["Initial_mM_Glucose"] == initial[0]) &
                          (data["Initial_mM_Acetate"] == initial[1])]
    fig, ax = plt.subplots()

    # Plot simulation results
    ax.plot(t, y[:, 0] * EXPERIMENT_VOLUME, color="b", label="Biomass")
    ax2 = plt.twinx(ax)
    ax2.plot(t, y[:, 1], color='r', label=f"Glucose (mM)")
    ax2.plot(t, y[:, 2], color='orange', label=f"Acetate (mM)")

    ax.set_ylabel('Biomass (g)', color='b')
    ax2.set_ylabel(f"Substrate (mM)", color='r')

    # Plot biomass from data
    mass_data = condition_data[condition_data["Type"] == "counts"]
    mass_data["Mass (g)"] = mass_data["Value"] * \
        MASS_PER_CELL * 1e-12  # pg to g
    mass_mean = mass_data.groupby("Time (h)")["Mass (g)"].mean().reset_index()
    mass_min = mass_data.groupby("Time (h)")["Mass (g)"].min().reset_index()
    mass_max = mass_data.groupby("Time (h)")["Mass (g)"].max().reset_index()

    ax.plot(mass_mean["Time (h)"], mass_mean["Mass (g)"], "b--")
    ax.fill_between(mass_min["Time (h)"], mass_min["Mass (g)"],
                    mass_max["Mass (g)"], color="b", alpha=0.2)

    # Plot substrates from data
    substrate_data = condition_data[condition_data["Type"]
                                    == "drawdown (umol)"]
    
    if substrate_data.size > 0:
        substrate_data["Drawdown (mM)"] = (substrate_data["Value"] *
                                        1e-3 /
                                        EXPERIMENT_VOLUME /
                                        substrate_data["Metabolite"].apply(lambda x: {"glucose": C_PER_GLUCOSE,
                                                                                        "acetate": C_PER_ACETATE}[x])
                                        )  # umol C to mM
        glucose_data = substrate_data[substrate_data["Metabolite"] == "glucose"]
        acetate_data = substrate_data[substrate_data["Metabolite"] == "acetate"]

        glucose_mean = glucose_data.groupby(
            "Time (h)")["Drawdown (mM)"].mean().reset_index()
        glucose_min = glucose_data.groupby(
            "Time (h)")["Drawdown (mM)"].min().reset_index()
        glucose_max = glucose_data.groupby(
            "Time (h)")["Drawdown (mM)"].max().reset_index()

        acetate_mean = acetate_data.groupby(
            "Time (h)")["Drawdown (mM)"].mean().reset_index()
        acetate_min = acetate_data.groupby(
            "Time (h)")["Drawdown (mM)"].min().reset_index()
        acetate_max = acetate_data.groupby(
            "Time (h)")["Drawdown (mM)"].max().reset_index()

        ax2.plot(glucose_mean["Time (h)"], glucose_mean["Drawdown (mM)"], "r--")
        ax2.fill_between(glucose_min["Time (h)"], glucose_min["Drawdown (mM)"],
                        glucose_max["Drawdown (mM)"], color="r", alpha=0.2)

        ax2.plot(acetate_mean["Time (h)"],
                acetate_mean["Drawdown (mM)"], "orange", linestyle="--")
        ax2.fill_between(acetate_min["Time (h)"], acetate_min["Drawdown (mM)"],
                        acetate_max["Drawdown (mM)"], color="orange", alpha=0.2)

    return fig, [ax, ax2]


def main():
    MODEL = "model/Rpom_05.xml"
    OUTDIR = "out/dFBA/CUE"
    BIOMASS_ID = "RPOM_provisional_biomass"

    DATA_FILE = "data/CUE/cue_data.csv"

    # Load data
    data = pd.read_csv(DATA_FILE)

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    # Load and set up model
    model = read_sbml_model(MODEL)
    setup_drawdown(model)

    # Set up Michaelis-Menten medium
    ex_glc = model.reactions.get_by_id("EX_glc")
    ex_ace = get_or_create_exchange(model, "ACET[e]")

    # Get V_maxes
    V_max_glc = abs(float(ex_glc._annotation["Experimental rate"]))
    V_max_ace = abs(float(ex_ace._annotation["Experimental rate"]))

    dynamic_medium = {
        ex_glc: MichaelisMentenBounds("Glucose[e]", V_max_glc, K_M),
        ex_ace: MichaelisMentenBounds("ACET[e]", V_max_ace, K_M)
    }

    # Turn off maintenance for now
    # TODO: bring back?
    # atpm = model.reactions.get_by_id("ATPM")
    # atpm.bounds = (0, 0)

    # Initial state
    initial_conditions = data[["Initial_mM_Glucose",
                               "Initial_mM_Acetate"]].drop_duplicates().values
    for initial_glucose, initial_acetate in initial_conditions:
        initial_biomass = (data[(data["Type"] == "counts") &
                                (data["Time (h)"] == 0) &
                                (data["Initial_mM_Glucose"] == initial_glucose) &
                                (data["Initial_mM_Acetate"] == initial_acetate)
                                ]["Value"].mean() * MASS_PER_CELL * 1e-12)
        tmax = (data[(data["Type"] == "counts") & (data["Initial_mM_Glucose"] == initial_glucose) & (
            data["Initial_mM_Acetate"] == initial_acetate)]["Time (h)"].max())

        initial = np.array([
            initial_biomass / EXPERIMENT_VOLUME,
            initial_glucose,
            initial_acetate
        ])

        t, y, l = dFBA(model,
                       BIOMASS_ID,
                       ["Glucose[e]", "ACET[e]"],
                       dynamic_medium,
                       initial,
                       tmax)

        # Plot data
        fig, _ = plot_result(t, y, [initial_glucose, initial_acetate],
                             V_max_ace, tmax, data)
        fig.set_size_inches(5, 3)
        fig.tight_layout()
        fig.savefig(os.path.join(
            OUTDIR, f"{initial_glucose:.2f}mM_glucose_{initial_acetate:.2f}mM_acetate_dFBA.png"))


if __name__ == "__main__":
    main()
