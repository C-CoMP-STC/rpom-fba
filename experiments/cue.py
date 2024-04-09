from matplotlib import pyplot as plt
import os
import pickle
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib
from cobra.io import read_sbml_model

from experiments.fast_dFBA import dFBA, make_shadow_price_listener, make_bge_listener, setup_drawdown, MichaelisMentenBounds, ConstantBounds, BoundFromData, plot_shadow_prices

from utils.cobra_utils import get_or_create_exchange
from utils.units import u
from model_building.model_factory import rebuild_and_get_model, DEFAULT_MODEL
from parameters.drawdown import *

matplotlib.use("Agg")

C_PER_GLUCOSE = 6
C_PER_ACETATE = 2


def plot_result(t, y, c, initial, condition_data, mass_units=True):
    glc_mw = 180.15588 * u.g / u.mol
    ace_mw = 59.04402 * u.g / u.mol

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

    # Plot bge
    ax3 = ax.twinx()
    ax3.plot(t, c, "k--")
    ax3.set_ylabel("BGE")
    ax3.set_ylim(0, 1)
    ax3.spines['right'].set_position(('outward', 50))

    ax.set_ylabel('Biomass (g/L)', color='b')
    ax2.set_ylabel(f"Substrate ({'g/L' if mass_units else 'mM'})", color='r')

    # Plot biomass from data
    mass_t = condition_data["mean"]["b_t"]
    mass_mean = condition_data["mean"]["b_s"]
    mass_min = condition_data["raw"]["raw_b"].min(axis=0)
    mass_max = condition_data["raw"]["raw_b"].max(axis=0)

    ax.plot(mass_t, mass_mean, "b--")
    ax.fill_between(mass_t, mass_min,
                    mass_max, color="b", alpha=0.2)
    
    # Plot substrates from data
    glucose_t = condition_data["raw"]["raw_g_t"]
    glucose = condition_data["raw"]["raw_g_s"]
    acetate_t = condition_data["raw"]["raw_a_t"]
    acetate = condition_data["raw"]["raw_a_s"]
    if mass_units:
        glucose = (glucose * u.mM * glc_mw).to("g/L").magnitude
        acetate = (acetate * u.mM * ace_mw).to("g/L").magnitude
    
    if glucose.size > 0:
        ax2.plot(glucose_t, glucose.mean(axis=0), "r--")
        ax2.fill_between(glucose_t, glucose.min(axis=0), glucose.max(axis=0), color="r", alpha=0.2)

    if acetate.size > 0:
        ax2.plot(acetate_t, acetate.mean(axis=0), "orange", linestyle="--")
        ax2.fill_between(acetate_t, acetate.min(axis=0), acetate.max(axis=0),
                        color="orange", alpha=0.2)

    ax2.legend()
    return fig, [ax, ax2]


def main():
    MODEL_BLUEPRINT = "model_building/blueprints/Rpom_05_ecoli_core.json"
    OUTDIR = "out/dFBA/CUE"
    # BIOMASS_ID = "RPOM_provisional_biomass"
    BIOMASS_ID = "Rpom_hwa_biomass"

    DATA_FILE = "data/clean/CUE/dFBA.pkl"

    # Load data
    with open(DATA_FILE, "rb") as f:
        data = pickle.load(f)

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)
    data_out = os.path.join(
        OUTDIR, datetime.now().strftime('%d-%b-%Y_%H:%M:%S'))
    os.makedirs(data_out)

    # Load and set up model
    # rebuild_and_get_model()  # MODEL_BLUEPRINT)
    # model = read_sbml_model("model/Rpom_05.xml")
    model = read_sbml_model("model/Rpom_05_hwa.xml")
    setup_drawdown(model)

    # Set up Michaelis-Menten medium
    ex_glc = model.reactions.get_by_id("EX_glc")
    ex_ace = get_or_create_exchange(model, "ACET[e]")

    # Get V_maxes
    V_max_glc = 10 * abs(float(ex_glc._annotation["Experimental rate"]))
    V_max_ace = 10 * abs(float(ex_ace._annotation["Experimental rate"]))

    for (initial_glucose, initial_acetate), condition_data in data.items():
        # Get biomass timeseries from data
        b_t, b = condition_data["mean"]["b_t"], condition_data["mean"]["b_s"]
        initial_biomass = b[0]
        tmax = b_t.max()

        initial = np.array([
            initial_biomass,
            initial_glucose.magnitude,
            initial_acetate.magnitude
        ])

        # Get glucose and acetate traces from data
        g_t = condition_data["mean"]["g_t"]
        g_s = condition_data["mean"]["g_s"]
        a_t = condition_data["mean"]["a_t"]
        a_s = condition_data["mean"]["a_s"]

        if initial_glucose == 0:
            g_t, g_s = np.array([0, tmax]), np.array([0, 0])
        if initial_acetate == 0:
            a_t, a_s = np.array([0, tmax]), np.array([0, 0])

        # Make bounds from data
        dynamic_medium = {
            ex_glc: BoundFromData("Glucose[e]", g_t, g_s, b_t, b),
            ex_ace: BoundFromData("ACET[e]", a_t, a_s, b_t, b),
        }

        t, y, l = dFBA(model,
                       BIOMASS_ID,
                       ["Glucose[e]", "ACET[e]"],
                       dynamic_medium,
                       CUE_VOLUME,
                       initial,
                       tmax,
                       terminate_on_infeasible=False,
                       listeners=[
                           make_shadow_price_listener(
                               model, ["Glucose[e]", "ACET[e]"], dynamic_medium),
                           make_bge_listener(
                               model, ["Glucose[e]", "ACET[e]"], dynamic_medium)
                       ],
                       dt=0.1)

        # Save data
        save_data = pd.DataFrame({"time (h)": t,
                                  "Biomass (g/L)": y[:, 0],
                                  "Glucose (mM)": y[:, 1],
                                  "Acetate (mM)": y[:, 2]})
        sp = pd.concat(l[0], axis=1, ignore_index=True).T
        save_data = pd.concat((save_data, sp), axis=1)
        save_data.to_csv(os.path.join(
            data_out, f"data_{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate.csv"), index=False)

        # Plot shadow prices
        fig, axs = plot_shadow_prices(l[0], t)
        axs[0].set_title(
            f"{initial_glucose.magnitude:.2f}mM Glucose, {initial_acetate.magnitude:.2f} mM Acetate")
        fig.set_size_inches(5, len(axs))
        fig.savefig(os.path.join(
            OUTDIR, f"{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_shadow_prices.png"))

        # Plot data
        c = l[1]
        fig, _ = plot_result(
            t, y, c, [initial_glucose, initial_acetate], condition_data, mass_units=False)
        fig.set_size_inches(6, 3)
        fig.tight_layout()
        fig.savefig(os.path.join(
            OUTDIR, f"{initial_glucose.magnitude:.2f}mM_glucose_{initial_acetate.magnitude:.2f}mM_acetate_dFBA.png"))


if __name__ == "__main__":
    main()
