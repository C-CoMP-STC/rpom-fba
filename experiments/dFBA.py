import json
import os

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from cobra.exceptions import Infeasible
from scipy.integrate import solve_ivp
from tqdm import tqdm
from parameters.fit_uptake_rates import get_mass_interpolator, michaelis_menten_dynamic_system

from utils.cobra_utils import get_or_create_exchange, set_active_bound, get_active_bound
from parameters.drawdown import *

MODEL = "model/Rpom_05.xml"
OUTDIR = "out/dFBA/"
CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"
UPTAKE_DATA_FILE = "data/drawdown_clean.csv"
GROWTH_DATA_FILE = "data/growth_curves_clean.csv"


def michaelis_menten_bounds(y, V_max, K_M):
    _, carbon_source = y  # expand the boundary species
    return V_max * carbon_source / (K_M + carbon_source)


def add_dynamic_bounds(y, exchange_rxn, V_max, K_M):
    """Use external concentrations to bound the uptake flux of glucose."""
    c_max_import = abs(michaelis_menten_bounds(y, V_max, K_M))

    set_active_bound(exchange_rxn, c_max_import)


def make_dynamic_system(model, biomass_rxn, exchange_rxn, K_M, pbar=None):
    exchange = model.reactions.get_by_id(exchange_rxn)
    V_max = get_active_bound(exchange)

    def dynamic_system(t, y):
        """Calculate the time derivative of external species."""

        biomass, _ = y  # expand the boundary species

        # Calculate the specific exchanges fluxes at the given external concentrations.
        with model:
            add_dynamic_bounds(y, exchange, V_max, K_M)

            # Ensure feasibility by setting feasibility as an objective,
            # then holding it fixed as a constraint
            cobra.util.add_lp_feasibility(model)
            cobra.util.fix_objective_as_constraint(model)

            # Using lexicographic optimization,
            # first optimize for biomass, then for the exchange flux
            # (holding optimal biomass as a constraint),
            # thus guaranteeing a unique optimal exchange flux.
            lex_constraints = cobra.util.add_lexicographic_constraints(
                model, [biomass_rxn, exchange_rxn], ['max', 'max'])
            fluxes = lex_constraints.values

        # Since the calculated fluxes are specific rates, we multiply them by the
        # biomass concentration to get the bulk exchange rates.
        fluxes[0] *= biomass
        # Fitted Vmax is specific, in units mM/hr/g, so need to multiply by g (not g/L)
        # to get mM/hr
        fluxes[1] *= biomass * COLONY_VOLUME * 1e-6

        # Update progress bar
        if dynamic_system.pbar is not None:
            dynamic_system.pbar.update(1)
            dynamic_system.pbar.set_description(
                f'{exchange_rxn} (t = {t:.3f})')

        return fluxes

    dynamic_system.pbar = pbar
    return dynamic_system


def make_infeasible_event(model, exchange_rxn, K_M, epsilon=1E-6):
    def infeasible_event(t, y):
        """
        Determine solution feasibility.

        Avoiding infeasible solutions is handled by solve_ivp's built-in event detection.
        This function re-solves the LP to determine whether or not the solution is feasible
        (and if not, how far it is from feasibility). When the sign of this function changes
        from -epsilon to positive, we know the solution is no longer feasible.

        """

        with model:
            exchange = model.reactions.get_by_id(exchange_rxn)
            V_max = get_active_bound(exchange)

            add_dynamic_bounds(y, exchange, V_max, K_M)

            cobra.util.add_lp_feasibility(model)
            feasibility = cobra.util.fix_objective_as_constraint(model)

        # solve_ivp event detection works by finding zeros, i.e. changes in sign.
        # however, if the problem is infeasible from the start, there will not
        # be a sign change.
        if t == 0:
            feasibility = 0

        return feasibility - epsilon

    return infeasible_event


def run_dFBA(model, exchange_rxn, y0, tmin, tmax, terminate_on_infeasible=True):
    # generate infeasibility detector
    infeasible_event = make_infeasible_event(
        model, exchange_rxn, K_M)
    infeasible_event.terminal = terminate_on_infeasible

    # TODO: solve_ivp supports a `vectorized` : Bool option,
    # indicating whether fun can be called in a vectorized fashion.
    #
    # so, make fun into a vectorized function if possible!

    with tqdm() as pbar:
        sol = solve_ivp(
            fun=make_dynamic_system(
                model, "RPOM_provisional_biomass", exchange_rxn, K_M, pbar),
            events=[infeasible_event],
            t_span=(tmin, tmax),
            y0=y0,
            rtol=1e-6,
            atol=1e-8,
            method='RK45'
        )

    return sol


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


def main():
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

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    for carbon_source, carbon_source_id in carbon_sources.items():
        if carbon_source_id != "Glucose[e]":
            continue

        with model:
            # Get Metabolite object and exchange reaction for the given carbon source
            exchange_rxn = get_or_create_exchange(model, carbon_source_id)

            vmax = abs(float(exchange_rxn._annotation["Experimental rate"]))

            set_active_bound(exchange_rxn, vmax)

            # Get initial conditions
            # TODO: set initial biomass from data?
            initial = uptake_data[uptake_data["Compound"] ==
                                  carbon_source]["InitialMetabolite_mM"].values[0]
            y0 = [COLONY_DRY_WEIGHT / COLONY_VOLUME,  # convert to concentration (g/L)
                  initial]

            # Get timespan of experiment
            tmin, tmax = 0, uptake_data[uptake_data["Compound"] ==
                                        carbon_source]["dt_hr"].values[0]

            data = run_dFBA(model, exchange_rxn.id, y0,
                            tmin, tmax)

            # Plot output
            fig, ax = plt.subplots()

            # Plot data
            ax.plot(data.t,
                    data.y.T[:, 0] * COLONY_VOLUME,
                    label="Biomass")
            ax2 = plt.twinx(ax)
            ax2.plot(data.t,
                     data.y.T[:, 1], color='r', label=f"[{carbon_source}] (mM)")

            ax.set_ylabel('Biomass (ug)', color='b')
            ax2.set_ylabel(carbon_source, color='r')

            # TODO: Plot fitted curve, points (?)
            col = f"{carbon_source}_predicted_mass"
            growth_on_carbon_source = growth_data[[col, "time (h)"]]
            growth_on_carbon_source = growth_on_carbon_source[~np.isnan(
                growth_on_carbon_source[col])]

            ax.plot(growth_on_carbon_source["time (h)"],
                    growth_on_carbon_source[col],
                    "b--",
                    label="Biomass (data)")

            # TODO: replace with fitted drawdown as in uptake rate fitting
            mass_curve = get_mass_interpolator(carbon_source, growth_data)
            t, x = michaelis_menten_dynamic_system(
                initial, mass_curve, -abs(vmax), K_M, tmax, dt=0.01)

            ax2.plot(t, x[:, 0], "r--",
                     label=f"Fitted {carbon_source_id} drawdown")

            # Save figure
            fig.set_size_inches(8, 6)
            fig.tight_layout()
            fig.savefig(os.path.join(OUTDIR, f"{carbon_source} dFBA.png"))


if __name__ == "__main__":
    main()
