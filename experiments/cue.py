from matplotlib import pyplot as plt
import os
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib

from experiments.fast_dFBA import dFBA, make_shadow_price_listener, setup_drawdown, MichaelisMentenBounds, ConstantBounds, plot_shadow_prices

from utils.cobra_utils import get_or_create_exchange
from utils.units import u
from model_building.model_factory import rebuild_and_get_model
from parameters.drawdown import *

matplotlib.use("Agg")

C_PER_GLUCOSE = 6
C_PER_ACETATE = 2


def plot_result(t, y, initial, data, mass_units=True):
    glc_mw = 180.15588 * u.g / u.mol
    ace_mw = 59.04402 * u.g / u.mol

    condition_data = data[(data["Initial_mM_Glucose"] == initial[0].magnitude) &
                          (data["Initial_mM_Acetate"] == initial[1].magnitude)]
    fig, ax = plt.subplots()

    # Plot simulation results
    ax.plot(t,
            y[:, 0],
            color="b",
            label="Biomass")
    ax2 = ax if mass_units else plt.twinx(ax)

    if mass_units:
        ax2.plot(t, ((y[:, 1] * u.mM) * glc_mw).to("g/L").magnitude,
                 color='r', label=f"Glucose (simulation)")
        ax2.plot(t, ((y[:, 2] * u.mM) * ace_mw).to("g/L").magnitude,
                 color='orange', label=f"Acetate (simulation)")
    else:
        ax2.plot(t, y[:, 1], color='r', label=f"Glucose (simulation)")
        ax2.plot(t, y[:, 2], color='orange', label=f"Acetate (simulation)")

    ax.set_ylabel('Biomass (g/L)', color='b')
    ax2.set_ylabel(f"Substrate ({'g/L' if mass_units else 'mM'})", color='r')

    # Plot biomass from data
    mass_data = condition_data[condition_data["Type"] == "counts"].copy()
    mass_data["Mass (g/L)"] = ((mass_data["Value"].values /
                                u.mL) * MASS_PER_CELL).to("g/L")
    mass_mean = mass_data.groupby(
        "Time (h)")["Mass (g/L)"].mean().reset_index()
    mass_min = mass_data.groupby("Time (h)")["Mass (g/L)"].min().reset_index()
    mass_max = mass_data.groupby("Time (h)")["Mass (g/L)"].max().reset_index()

    ax.plot(mass_mean["Time (h)"], mass_mean["Mass (g/L)"], "b--")
    ax.fill_between(mass_min["Time (h)"], mass_min["Mass (g/L)"],
                    mass_max["Mass (g/L)"], color="b", alpha=0.2)

    # Plot substrates from data
    substrate_data = condition_data[condition_data["Type"]
                                    == "drawdown (umol)"].copy()

    if substrate_data.size > 0:
        substrate_data["Drawdown (mM)"] = ((substrate_data["Value"].values * u.umol).to("mmol") /
                                           CUE_VOLUME /
                                           substrate_data["Metabolite"].apply(lambda x: {"glucose": C_PER_GLUCOSE,
                                                                                         "acetate": C_PER_ACETATE}[x]).values
                                           ).to("mM").magnitude
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

        if mass_units:
            glucose_mean["Drawdown (mM)"] = (
                (glucose_mean["Drawdown (mM)"].values * u.mM) * glc_mw).to("g/L").magnitude
            glucose_min["Drawdown (mM)"] = (
                (glucose_min["Drawdown (mM)"].values * u.mM) * glc_mw).to("g/L").magnitude
            glucose_max["Drawdown (mM)"] = (
                (glucose_max["Drawdown (mM)"].values * u.mM) * glc_mw).to("g/L").magnitude

            acetate_mean["Drawdown (mM)"] = (
                (acetate_mean["Drawdown (mM)"].values * u.mM) * ace_mw).to("g/L").magnitude
            acetate_min["Drawdown (mM)"] = (
                (acetate_min["Drawdown (mM)"].values * u.mM) * ace_mw).to("g/L").magnitude
            acetate_max["Drawdown (mM)"] = (
                (acetate_max["Drawdown (mM)"].values * u.mM) * ace_mw).to("g/L").magnitude

        ax2.plot(glucose_mean["Time (h)"],
                 glucose_mean["Drawdown (mM)"], "r--")
        ax2.fill_between(glucose_min["Time (h)"], glucose_min["Drawdown (mM)"],
                         glucose_max["Drawdown (mM)"], color="r", alpha=0.2)

        ax2.plot(acetate_mean["Time (h)"],
                 acetate_mean["Drawdown (mM)"], "orange", linestyle="--")
        ax2.fill_between(acetate_min["Time (h)"], acetate_min["Drawdown (mM)"],
                         acetate_max["Drawdown (mM)"], color="orange", alpha=0.2)

    ax2.legend()

    return fig, [ax, ax2]


def main():
    MODEL_BLUEPRINT = "model_building/blueprints/Rpom_05_ecoli_core.json"
    OUTDIR = "out/dFBA/CUE"
    BIOMASS_ID = "RPOM_provisional_biomass"

    DATA_FILE = "data/clean/CUE/cue_data.csv"

    # Load data
    data = pd.read_csv(DATA_FILE)

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)
    data_out = os.path.join(
        OUTDIR, datetime.now().strftime('%d-%b-%Y_%H:%M:%S'))
    os.makedirs(data_out)

    # Load and set up model
    model = rebuild_and_get_model()# MODEL_BLUEPRINT)
    setup_drawdown(model)

    # Set up Michaelis-Menten medium
    ex_glc = model.reactions.get_by_id("EX_glc")
    ex_ace = get_or_create_exchange(model, "ACET[e]")

    # Get V_maxes
    V_max_glc = 10 * abs(float(ex_glc._annotation["Experimental rate"]))
    V_max_ace = 10 * abs(float(ex_ace._annotation["Experimental rate"]))

    dynamic_medium = {
        ex_glc: MichaelisMentenBounds("Glucose[e]", V_max_glc, K_M.to("mM").magnitude),
        ex_ace: MichaelisMentenBounds(
            "ACET[e]", V_max_ace, K_M.to("mM").magnitude)
    }

    # Initial state
    initial_conditions = data[["Initial_mM_Glucose",
                               "Initial_mM_Acetate"]].drop_duplicates().values
    for initial_glucose, initial_acetate in initial_conditions:
        
        initial_biomass = (data[(data["Type"] == "counts") &
                                (data["Time (h)"] == 0) &
                                (data["Initial_mM_Glucose"] == initial_glucose) &
                                (data["Initial_mM_Acetate"] == initial_acetate)
                                ]["Value"].mean() * (1/u.mL) * MASS_PER_CELL).to("g/L")
        tmax = (data[(data["Type"] == "counts") & (data["Initial_mM_Glucose"] == initial_glucose) & (
            data["Initial_mM_Acetate"] == initial_acetate)]["Time (h)"].max())

        initial_glucose *= u.mM
        initial_acetate *= u.mM

        initial = np.array([
            initial_biomass.magnitude,
            initial_glucose.magnitude,
            initial_acetate.magnitude
        ])

        t, y, l = dFBA(model,
                       BIOMASS_ID,
                       ["Glucose[e]", "ACET[e]"],
                       dynamic_medium,
                       CUE_VOLUME,
                       initial,
                       tmax,
                       listeners=[make_shadow_price_listener(
                           model, ["Glucose[e]", "ACET[e]"], dynamic_medium)],
                       dt=0.1)

        # Save data
        save_data = pd.DataFrame({"time (h)": t,
                                  "Biomass (g/L)": y[:, 0],
                                  "Glucose (mM)": y[:, 1],
                                  "Acetate (mM)": y[:, 2]})
        sp = pd.concat(l, axis=1, ignore_index=True).T
        save_data = pd.concat((save_data, sp), axis=1)
        save_data.to_csv(os.path.join(
            data_out, f"data_{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate.csv"), index=False)

        # Plot data
        fig, axs = plot_shadow_prices(l, t)
        axs[0].set_title(
            f"{initial_glucose.magnitude:.2f}mM Glucose, {initial_acetate.magnitude:.2f} mM Acetate")
        fig.set_size_inches(5, len(axs))
        fig.savefig(os.path.join(
            OUTDIR, f"{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_shadow_prices.png"))

        fig, _ = plot_result(t, y, [initial_glucose, initial_acetate], data, mass_units=False)
        fig.set_size_inches(5, 3)
        fig.tight_layout()
        fig.savefig(os.path.join(
            OUTDIR, f"{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_dFBA.png"))


if __name__ == "__main__":
    main()
