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


EXPERIMENT_VOLUME = 100  # mL


def plot_result(t, y, initial, V_max, t_max, data):
    fig, ax = plt.subplots()

    # Plot data
    ax.plot(t, y[:, 0] * EXPERIMENT_VOLUME, color="b", label="Biomass")
    ax2 = plt.twinx(ax)
    ax2.plot(t, y[:, 1], color='r', label=f"Glucose (mM)")
    ax2.plot(t, y[:, 2], color='orange', label=f"Acetate (mM)")

    ax.set_ylabel('Biomass (g)', color='b')
    ax2.set_ylabel(f"Substrate (mM)", color='r')

    # col = f"{carbon_source}_predicted_mass"
    # growth_on_carbon_source = growth_data[[col, "time (h)"]]
    # growth_on_carbon_source = growth_on_carbon_source[~np.isnan(
    #     growth_on_carbon_source[col])]

    # ax.plot(growth_on_carbon_source["time (h)"],
    #         growth_on_carbon_source[col] * 1e6,  # Convert g -> ug
    #         "b--",
    #         label="Biomass (data)")

    # mass_curve = get_mass_interpolator(carbon_source, growth_data)
    # t, x = michaelis_menten_dynamic_system(
    #     initial_C, mass_curve, -abs(V_max), K_M, t_max, dt=0.01)

    # ax2.plot(t, x[:, 0], "r--",
    #          label=f"Fitted {carbon_source} drawdown")

    return fig, [ax, ax2]


def main():
    MODEL = "model/Rpom_05.xml"
    OUTDIR = "out/dFBA/CUE"
    BIOMASS_ID = "RPOM_provisional_biomass"

    DATA_FILE = "data/CUE/cue_data.csv"

    # TODO: Estimated from E coli, get a better number if possible and move to a single source
    MASS_PER_CELL = 0.95  # pg

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
