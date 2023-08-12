from itertools import pairwise
from utils.cobra_utils import get_or_create_exchange, set_active_bound
from parameters.fit_uptake_rates import (get_mass_interpolator,
                                         michaelis_menten_dynamic_system)
from parameters.drawdown import *
import matplotlib.pyplot as plt
import json
import os

import cobra
import matplotlib
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from tqdm import tqdm

matplotlib.use("Agg")


def rk45(df_dt, y0, tmin, tmax, dt=0.01, terminate_on_error=True, pbar=True, pbar_desc=None, listeners=None):
    def rk45_step(df_dt, y0, dt):
        k1 = df_dt(y0) * dt
        k2 = df_dt(y0 + 0.5 * k1) * dt
        k3 = df_dt(y0 + 0.5 * k2) * dt
        k4 = df_dt(y0 + k3) * dt

        return y0 + (k1 + 2*k2 + 2*k3 + k4)/6

    t_range = np.arange(tmin, tmax, dt)
    result = np.zeros((t_range.size, y0.size))
    result[0, :] = y0
    listener_data = ([listener(y0) for listener in listeners]
                     if listeners is not None else [])

    t_index = range(1, len(t_range))
    for i in tqdm(t_index, pbar_desc) if pbar else t_index:
        try:
            y = rk45_step(df_dt, result[i-1], dt)
            result[i, :] = y

            # Run listeners
            listener_data += ([listener(y) for listener in listeners]
                              if listeners is not None else [])

        except Exception as e:
            if terminate_on_error:
                return t_range[:i], result[:i, :], listener_data
            raise e

    return t_range, result, listener_data


def dFBA(model, biomass_id, dynamic_medium, V_max, K_M, y0, tmax, dt=0.01, terminate_on_infeasible=True, listeners=None, desc=""):
    medium_ids = [rxn.id for rxn in dynamic_medium.keys()]

    def df_dt(y):
        biomass, carbon_source = y

        with model:
            for exchange, updater in dynamic_medium.items():
                if updater == "michaelis-menten":
                    mm_bound = abs(V_max * carbon_source / (K_M + carbon_source))
                    set_active_bound(exchange, mm_bound)

            # Using lexicographic optimization,
            # first optimize for biomass, then for the exchange fluxes
            # (holding optimal biomass as a constraint),
            # thus guaranteeing a unique optimal set of exchange fluxes.
            lex_constraints = cobra.util.add_lexicographic_constraints(
                model, [biomass_id] + medium_ids, ['max' for _ in range(len(medium_ids) + 1)])
            fluxes = lex_constraints.values

        # Fluxes are specific rates, so we multiply them by the
        # biomass concentration to get the bulk exchange rates.
        # Fitted Vmax is also specific (mM/hr/g), so we further multiply
        # by volume to get final units of mM/hr
        fluxes *= biomass
        fluxes[1] *= COLONY_VOLUME * 1e-6

        return fluxes

    return rk45(df_dt, y0, 0, tmax, dt, terminate_on_infeasible, pbar_desc=desc, listeners=listeners)


def make_shadow_price_listener(model, V_max, C_exchange_id, n=10):
    exchange = model.reactions.get_by_id(C_exchange_id)

    def shadow_price_listener(y):
        _, carbon_source = y

        with model:
            mm_bound = abs(V_max * carbon_source / (K_M + carbon_source))
            set_active_bound(exchange, mm_bound)

            sol = model.optimize()
            max_shadow_price_metabolites = sol.shadow_prices.abs(
            ).sort_values(ascending=False)[:n].index

            return sol.shadow_prices[max_shadow_price_metabolites]

    return shadow_price_listener


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
            initial_biomass = growth_data[f"{carbon_source}_predicted_mass"].values[0] / (
                COLONY_VOLUME * 1e-6)
            y0 = np.array([initial_biomass,
                           initial_C])

            # Run experiment
            tmax = uptake_data[uptake_data["Compound"] ==
                               carbon_source]["dt_hr"].values[0]

            dynamic_medium = {exchange_rxn : "michaelis-menten"}
            t, y, listeners = dFBA(model, BIOMASS_ID, dynamic_medium, V_max,
                                   K_M, y0, tmax, dt=0.01, terminate_on_infeasible=True,
                                   listeners=[make_shadow_price_listener(model, V_max, exchange_rxn.id)],
                                   desc=carbon_source)

            # Plot shadow prices over time
            all_metabolites = set()
            for shadow_prices in listeners:
                all_metabolites.update(shadow_prices.index)
            df = pd.DataFrame({"time": t})
            for metabolite in all_metabolites:
                df[metabolite] = [shadow_prices[metabolite]
                                  if metabolite in shadow_prices.index else 0 for shadow_prices in listeners]
            sum_abs = df.loc[:, df.columns != 'time'].abs().sum(axis=1)
            max_abs = df.loc[:, df.columns != 'time'].abs().max(axis=1)
            df["sum_abs"] = sum_abs
            df["max_abs"] = max_abs

            fig, ax = plt.subplots()
            bottom = np.zeros(df.shape[0])
            for i, metabolite in enumerate(all_metabolites):
                heights = (df[metabolite].abs() / df["sum_abs"]).values
                pc = matplotlib.collections.PatchCollection(
                    [
                        matplotlib.patches.Rectangle(
                            (df.iloc[r1]["time"], bottom[r1]),
                            df.iloc[r2]["time"] - df.iloc[r1]["time"],
                            heights[r1])
                        for r1, r2 in pairwise(range(df.shape[0]))
                    ]
                )
                pc.set_cmap("Spectral")
                pc.set_norm(plt.Normalize(-max(max_abs), max(max_abs)))
                pc.set_array(df[metabolite])
                ax.add_collection(pc)

                ax.plot(df["time"], bottom + heights,
                        "w", label=f"{i}: {metabolite}")

                max_height = heights.max()
                max_height_idx = heights.argmax()
                ax.text(df["time"].values[max_height_idx], bottom[max_height_idx],
                        f"{i}", horizontalalignment="left", verticalalignment="bottom")
                bottom += heights

            leg = fig.legend(handlelength=0, handletextpad=0,
                             loc="center left", bbox_to_anchor=(1, 0.5))
            for item in leg.legendHandles:
                item.set_visible(False)
            ax.set_xlim(t.min(), t.max())
            ax.set_xlabel("Time (hr)")
            ax.set_ylabel("Shadow Price (normalized magnitude)")
            ax.set_title(f"Shadow prices over time on {carbon_source}")
            fig.colorbar(pc, location="left", pad=0.2)
            fig.savefig(os.path.join(
                OUTDIR, f"{carbon_source} shadow prices.png"), bbox_inches='tight')

            # Plot data
            fig, _ = plot_data(t, y, carbon_source, initial_C,
                               V_max, tmax, growth_data)
            fig.set_size_inches(5, 3)
            fig.tight_layout()
            fig.savefig(os.path.join(OUTDIR, f"{carbon_source} dFBA.png"))


if __name__ == "__main__":
    main()
