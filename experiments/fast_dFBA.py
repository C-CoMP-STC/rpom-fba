from parameters.drawdown import *
from parameters.fit_uptake_rates import get_mass_interpolator, michaelis_menten_dynamic_system
from utils.cobra_utils import get_or_create_exchange, set_active_bound
from tqdm import tqdm
from cobra.io import read_sbml_model
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import os

import cobra
import matplotlib
matplotlib.use("Agg")


def rk45(df_dt, y0, tmin, tmax, dt=0.01, terminate_on_error=True, pbar=True, pbar_desc=None):
    def rk45_step(df_dt, y0, dt):
        k1 = df_dt(y0) * dt
        k2 = df_dt(y0 + 0.5 * k1) * dt
        k3 = df_dt(y0 + 0.5 * k2) * dt
        k4 = df_dt(y0 + k3) * dt

        return y0 + (k1 + 2*k2 + 2*k3 + k4)/6

    t_range = np.arange(tmin, tmax, dt)
    result = np.zeros((t_range.size, y0.size))
    result[0, :] = y0

    t_index = range(1, len(t_range))
    for i in tqdm(t_index, pbar_desc) if pbar else t_index:
        try:
            result[i, :] = rk45_step(df_dt, result[i-1], dt)
        except Exception as e:
            if terminate_on_error:
                return t_range[:i], result[:i, :]
            raise e

    return t_range, result


def dFBA(model, biomass_id, C_exchange_id, V_max, K_M, y0, tmax, dt=0.01, terminate_on_infeasible=True):
    exchange = model.reactions.get_by_id(C_exchange_id)

    def df_dt(y):
        biomass, carbon_source = y

        with model:
            mm_bound = abs(V_max * carbon_source / (K_M + carbon_source))
            set_active_bound(exchange, mm_bound)

            lex_constraints = cobra.util.add_lexicographic_constraints(
                model, [biomass_id, C_exchange_id], ['max', 'max'])
            fluxes = lex_constraints.values

        # Fluxes are specific rates, so we multiply them by the
        # biomass concentration to get the bulk exchange rates.
        # Fitted Vmax is also specific (mM/hr/g), so we further multiply
        # by volume to get final units of mM/hr
        fluxes *= biomass
        fluxes[1] *= COLONY_VOLUME * 1e-6

        return fluxes

    return rk45(df_dt, y0, 0, tmax, dt, terminate_on_infeasible, pbar_desc=C_exchange_id)


def setup_drawdown(model):
    # Growth is infeasible on the seawater medium as it currently is,
    # needs to be supplemented with FE+2 (also increase everything to 1000 to not be limiting)
    supp_medium = {k: 1000. for k, v in model.medium.items()}
    supp_medium["EX_fe2"] = 1000.
    model.medium = supp_medium

    # Remove biotin from objective temporarily as biotin is blocking
    # TODO: fix biotin production?
    biotin = model.metabolites.get_by_id("BIOTIN[c]")
    biomass = model.reactions.get_by_id("RPOM_provisional_biomass")
    biomass.subtract_metabolites({biotin: biomass.metabolites[biotin]})

    # TODO: Growth is currently slow on glucose -
    # try changing maintenance requirement
    # atpm = model.reactions.get_by_id("ATPM")
    # atpm.lower_bound = atpm.upper_bound = 10

    # TODO: Growth is currently slow on glucose -
    # testing lower peptidoglycan requirement
    # murein = model.metabolites.get_by_id("CPD0-2278[p]")
    # biomass.add_metabolites({murein: 1/10 * abs(biomass.metabolites[murein])})


def plot_data(t, y, carbon_source, initial_C, V_max, t_max, growth_data):
    fig, ax = plt.subplots()

    # Plot data
    ax.plot(t, y[:, 0] * COLONY_VOLUME, color="b", label="Biomass")
    ax2 = plt.twinx(ax)
    ax2.plot(t, y[:, 1], color='r', label=f"[{carbon_source}] (mM)")

    ax.set_ylabel('Biomass (ug)', color='b')
    ax2.set_ylabel(f"{carbon_source} (mM)", color='r')

    col = f"{carbon_source}_predicted_mass"
    growth_on_carbon_source = growth_data[[col, "time (h)"]]
    growth_on_carbon_source = growth_on_carbon_source[~np.isnan(
        growth_on_carbon_source[col])]

    ax.plot(growth_on_carbon_source["time (h)"],
            growth_on_carbon_source[col] * 1e6,  # Convert g -> ug
            "b--",
            label="Biomass (data)")

    mass_curve = get_mass_interpolator(carbon_source, growth_data)
    t, x = michaelis_menten_dynamic_system(
        initial_C, mass_curve, -abs(V_max), K_M, t_max, dt=0.01)

    ax2.plot(t, x[:, 0], "r--",
             label=f"Fitted {carbon_source} drawdown")

    return fig, [ax, ax2]


def main():
    MODEL = "model/Rpom_05.xml"
    BIOMASS_ID = "RPOM_provisional_biomass"
    CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"
    UPTAKE_DATA_FILE = "data/drawdown_clean.csv"
    GROWTH_DATA_FILE = "data/growth_curves_clean.csv"
    OUTDIR = "out/dFBA/"

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    model = read_sbml_model(MODEL)
    setup_drawdown(model)

    # Load carbon sources to test
    with open(CARBON_SOURCES, "r") as f:
        carbon_sources = json.load(f)
    carbon_sources = {
        k: v
        for k, v in carbon_sources.items()
        if len(model.metabolites.query(lambda m: m.id == v)) == 1
    }

    # Load data for compoarison
    uptake_data = pd.read_csv(UPTAKE_DATA_FILE)
    growth_data = pd.read_csv(GROWTH_DATA_FILE)

    for carbon_source, carbon_source_id in carbon_sources.items():
        with model:
            # Get Metabolite object and exchange reaction for the given carbon source
            exchange_rxn = get_or_create_exchange(model, carbon_source_id)
            V_max = abs(float(exchange_rxn._annotation["Experimental rate"]))

            # Get initial conditions
            # TODO: set initial biomass from data?
            initial_C = uptake_data[uptake_data["Compound"] ==
                                    carbon_source]["InitialMetabolite_mM"].values[0]
            initial_biomass = growth_data[f"{carbon_source}_predicted_mass"].values[0] / (COLONY_VOLUME * 1e-6)
            y0 = np.array([initial_biomass,
                           initial_C])

            # Run experiment
            tmax = uptake_data[uptake_data["Compound"] ==
                               carbon_source]["dt_hr"].values[0]

            t, y = dFBA(model, BIOMASS_ID, exchange_rxn.id, V_max,
                        K_M, y0, tmax, dt=0.01, terminate_on_infeasible=True)

            # Plot data
            fig, _ = plot_data(t, y, carbon_source, initial_C,
                               V_max, tmax, growth_data)
            fig.set_size_inches(6, 4)
            fig.tight_layout()
            fig.savefig(os.path.join(OUTDIR, f"{carbon_source} dFBA.png"))


if __name__ == "__main__":
    main()
