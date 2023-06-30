import json
import os

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from scipy.integrate import solve_ivp
from tqdm import tqdm

from utils.cobra_utils import change_compartment, get_or_create_exchange, set_active_bound

MODEL = "model/Rpom_05.xml"
OUTDIR = "out/dFBA/"
CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"
UPTAKE_RATES = "parameters/uptake_rates/fitted_uptake_rates.json"
UPTAKE_DATA_FILE = "data/drawdown_clean.csv"


COLONY_VOLUME = 220  # uL
COLONY_DRY_WEIGHT = 405.3053598  # ug


def michaelis_menten_bounds(y, V_max, K_M):
    _, carbon_source = y # expand the boundary species
    return V_max * carbon_source / (K_M + carbon_source)


def add_dynamic_bounds(y, exchange_rxn, V_max, K_M):
    """Use external concentrations to bound the uptake flux of glucose."""
    c_max_import = michaelis_menten_bounds(y, V_max, K_M)

    set_active_bound(exchange_rxn, c_max_import)


def make_dynamic_system(model, biomass_rxn, exchange_rxn, K_M, pbar=None):
    exchange = model.reactions.get_by_id(exchange_rxn)
    V_max = exchange_rxn.lower_bound

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

        # Since the calculated fluxes are specific rates, we multiply them by the
        # biomass concentration to get the bulk exchange rates.
        fluxes = lex_constraints.values
        fluxes *= biomass

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
            V_max = model.reactions.get_by_id(exchange_rxn).lower_bound

            add_dynamic_bounds(model, y, exchange_rxn, V_max, K_M)

            cobra.util.add_lp_feasibility(model)
            feasibility = cobra.util.fix_objective_as_constraint(model)

        # solve_ivp event detection works by finding zeros, i.e. changes in sign.
        # however, if the problem is infeasible from the start, there will not
        # be a sign change.
        if t == 0:
            feasibility = 0

        print(feasibility)

        return feasibility - epsilon

    return infeasible_event


def experiment(model, exchange_rxn, y0, tmin, tmax):
    # TODO: Extremely arbitrary!!!
    K_M = 5

    # generate infeasibility detector
    infeasible_event = make_infeasible_event(model, exchange_rxn, K_M)

    # Terminate at the first instance of infeasibility
    infeasible_event.terminal = True

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
            # t_eval=T,
            rtol=1e-6,
            atol=1e-8,
            method='RK45'
        )

    return sol


def plot_data(ivp_solution, carbon_source, ax):
    ax.plot(ivp_solution.t, ivp_solution.y.T[:, 0], label="Biomass")
    ax2 = plt.twinx(ax)
    ax2.plot(ivp_solution.t,
             ivp_solution.y.T[:, 1], color='r', label=carbon_source)

    ax.set_ylabel('Biomass', color='b')
    ax2.set_ylabel(carbon_source, color='r')

    return ax, ax2


def main():
    model = read_sbml_model(MODEL)

    # Growth is infeasible on the seawater medium as it currently is,
    # supplement with the following:
    # supplement = ['CO+2[c]', 'CU+2[c]', 'MN+2[c]', 'CPD-3[c]', 'NI+2[c]', 'THIAMINE[c]', "FE+2[c]", "ZN+2[c]"]
    # supplement.append("Pi[c]")
    # supp_medium = model.medium
    # for met in supplement:
    #     ex = get_or_create_exchange(model, met)
    #     supp_medium[ex.id] = 1000
    # model.medium = supp_medium
    supp_medium = {k : 1000. for k, v in model.medium.items()}
    supp_medium["EX_fe2"] = 1000.
    model.medium = supp_medium

    # Remove biotin from objective temporarily as biotin is blocking
    # TODO: fix biotin production?
    biotin = model.metabolites.get_by_id("BIOTIN[c]")
    biomass = model.reactions.get_by_id("RPOM_provisional_biomass")
    biomass.subtract_metabolites({biotin: biomass.metabolites[biotin]})

    # Load carbon sources to test
    with open(CARBON_SOURCES, "r") as f:
        carbon_sources = json.load(f)
    carbon_sources = {
        k: v
        for k, v in carbon_sources.items()
        if len(model.metabolites.query(lambda m: m.id == v)) == 1
    }

    # Load uptake rate source data for comparison
    uptake_data = pd.read_csv(UPTAKE_DATA_FILE)

    with open(UPTAKE_RATES, "r") as f:
        uptake_rates = json.load(f)

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    for carbon_source, carbon_source_id in carbon_sources.items():
        if carbon_source_id != "Glucose[e]":
            continue

        with model:
            # Get Metabolite object and exchange reaction for the given carbon source
            exchange_rxn = get_or_create_exchange(model, carbon_source_id)

            # TODO: Feeding extra for now, due to infeasibility on calculated
            # feeding rate (at least for glucose)
            rate = 2 * \
                abs(float(exchange_rxn._annotation["Experimental rate"]))

            set_active_bound(exchange_rxn, rate)

            # Run the experiment, get data
            # Biomass, carbon source
            initial = uptake_data[uptake_data["Compound"] ==
                                  carbon_source]["InitialMetabolite_mM"].values[0]
            # convert to mmol/gDCW
            initial *= COLONY_VOLUME / COLONY_DRY_WEIGHT

            # TODO: set initial biomass from data?
            y0 = [COLONY_DRY_WEIGHT * 1e-6, initial]

            tmin, tmax = 0, 24
            data = experiment(model, exchange_rxn.id, y0, tmin, tmax)

            # Plot output
            fig, ax = plt.subplots()

            # Plot depletion of carbon source over 24 hrs
            _, ax2 = plot_data(data, carbon_source_id, ax)

            # TODO: Plot fitted curve, points (?)

            rate = uptake_rates[carbon_source]
            t = np.linspace(tmin, tmax, 50)
            ax2.plot(t, initial * np.exp(rate * t),
                     "r--", label="Theoretical drawdown")

            fig.legend()

            # Save figure
            fig.set_size_inches(8, 6)
            fig.tight_layout()
            fig.savefig(os.path.join(OUTDIR, f"{carbon_source} dFBA.png"))


if __name__ == "__main__":
    main()
