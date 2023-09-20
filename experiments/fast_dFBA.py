import matplotlib.pyplot as plt
import json
import os
from itertools import pairwise

import cobra
import matplotlib
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from tqdm import tqdm

from parameters.drawdown import *
from parameters.fit_uptake_rates import michaelis_menten_dynamic_system
from utils.cobra_utils import get_or_create_exchange, set_active_bound
from utils.math import get_interpolator, rk45
from utils.units import u

matplotlib.use("Agg")


class MichaelisMentenBounds:
    def __init__(self, substrate_id, V_max, K_M):
        self.substrate_id = substrate_id
        self.V_max = V_max
        self.K_M = K_M

    def bound(self, exchange, concentration):
        concentration = max(concentration, 0 * u.mM)
        mm_bound = abs(self.V_max * concentration / (K_M + concentration))
        set_active_bound(exchange, mm_bound.magnitude)


def dFBA(model, biomass_id, substrate_ids, dynamic_medium, volume, y0, tmax, dt=0.01, terminate_on_infeasible=True, listeners=None, desc=""):
    medium_ids = [rxn.id for rxn in dynamic_medium.keys()]

    def df_dt(y):
        biomass = y[0] * u.g/u.L
        substrates = dict(zip(substrate_ids, y[1:] * u.mM))

        with model:
            for exchange, bounds in dynamic_medium.items():
                bounds.bound(exchange, substrates[bounds.substrate_id])

            # Using lexicographic optimization,
            # first optimize for biomass, then for the exchange fluxes
            # (holding optimal biomass as a constraint),
            # thus guaranteeing a unique optimal set of exchange fluxes.
            lex_constraints = cobra.util.add_lexicographic_constraints(
                model, [biomass_id] + medium_ids, ['max' for _ in range(len(medium_ids) + 1)])
            # TODO: FIX!!!! biomass is not in fx units
            biomass_flux = lex_constraints.values[0] * (1/u.hr)
            medium_fluxes = lex_constraints.values[1:] * u.fx

        # Fluxes are specific rates, so we multiply them by the
        # biomass concentration to get the bulk exchange rates.
        # Fitted Vmax is also specific (mM/hr/g), so we further multiply
        # by volume to get final units of mM/hr
        fluxes *= biomass
        fluxes = fluxes.to("mmol/L/hr").magnitude

        # TODO: WRONG!!! Pass in a volume
        fluxes[1:] *= volume.to("L").magnitude

        return fluxes

    return rk45(df_dt, y0, 0, tmax, dt, terminate_on_infeasible, pbar_desc=desc, listeners=listeners)


def make_shadow_price_listener(model, substrate_ids, dynamic_medium, n=10):
    def shadow_price_listener(y):
        substrates = dict(zip(substrate_ids, y[1:] * u.mM))

        with model:
            for exchange, bounds in dynamic_medium.items():
                bounds.bound(exchange, substrates[bounds.substrate_id])

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
    supp_medium["EX_o2"] = 65.
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
    ax.plot(t,
            (y[:, 0] * u.g/u.L * COLONY_VOLUME.to("L")).to("ug"),
            color="b",
            label="Biomass")
    ax2 = plt.twinx(ax)
    ax2.plot(t, y[:, 1], color='r', label=f"[{carbon_source}] (mM)")

    ax.set_ylabel('Biomass (ug)', color='b')
    ax2.set_ylabel(f"{carbon_source} (mM)", color='r')

    col = f"{carbon_source}_predicted_mass"
    growth_on_carbon_source = growth_data[[col, "time (h)"]]
    growth_on_carbon_source = growth_on_carbon_source[~np.isnan(
        growth_on_carbon_source[col])]

    ax.plot(growth_on_carbon_source["time (h)"],
            (growth_on_carbon_source[col].values * u.g).to("ug"),
            "b--",
            label="Biomass (data)")

    t = growth_data["time (h)"] * u.h
    y = growth_data[f"{carbon_source}_predicted_mass"] * u.g
    mass_curve = get_interpolator(t, y)
    t, x = michaelis_menten_dynamic_system(
        initial_C.magnitude, mass_curve, -abs(V_max), K_M.magnitude, t_max, dt=0.01)

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
            initial_C = uptake_data[uptake_data["Compound"] ==
                                    carbon_source]["InitialMetabolite_mM"].values[0] * u.mM
            initial_biomass = growth_data[f"{carbon_source}_predicted_mass"].values[0] * u.g / (
                COLONY_VOLUME.to("L"))
            y0 = np.array([initial_biomass.magnitude,
                           initial_C.magnitude])

            # Run experiment
            tmax = uptake_data[uptake_data["Compound"] ==
                               carbon_source]["dt_hr"].values[0]

            dynamic_medium = {exchange_rxn: MichaelisMentenBounds(
                carbon_source_id, V_max, K_M)}
            t, y, listeners = dFBA(model, BIOMASS_ID, [carbon_source_id], dynamic_medium, COLONY_VOLUME,
                                   y0, tmax, dt=0.01, terminate_on_infeasible=True,
                                   listeners=[make_shadow_price_listener(
                                       model, [carbon_source_id], dynamic_medium)],
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
