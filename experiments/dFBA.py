import json
import os

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from scipy.integrate import solve_ivp
from tqdm import tqdm

MODEL = "model/Rpom_05.xml"
OUTDIR = "out/dFBA/"
CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"
UPTAKE_RATES = "parameters/uptake_rates/fitted_uptake_rates.json"
UPTAKE_DATA_FILE = "data/drawdown_clean.csv"

# TODO: Let's output and save data for all exchange fluxes in addition to biomass and external C


def add_dynamic_bounds(model, y, exchange_rxn, V_max, K_M):
    """Use external concentrations to bound the uptake flux of glucose."""
    biomass, carbon_source = y  # expand the boundary species
    c_max_import = V_max * carbon_source / (K_M + carbon_source)
    if c_max_import > 0:
        c_max_import = 0
    model.reactions.get_by_id(exchange_rxn).lower_bound = c_max_import


def make_dynamic_system(model, biomass_rxn, exchange_rxn, K_M, pbar=None):
    V_max = model.reactions.get_by_id(exchange_rxn).lower_bound

    def dynamic_system(t, y):
        """Calculate the time derivative of external species."""

        biomass, carbon_source = y  # expand the boundary species

        # Calculate the specific exchanges fluxes at the given external concentrations.
        with model:
            add_dynamic_bounds(model, y, exchange_rxn, V_max, K_M)

            # Check feasbility by adding a feasibility constraint(?)
            cobra.util.add_lp_feasibility(model)
            feasibility = cobra.util.fix_objective_as_constraint(model)

            # Optimize biomass, and then carbon source exchange (why?)
            lex_constraints = cobra.util.add_lexicographic_constraints(
                model, [biomass_rxn, exchange_rxn], ['max', 'max'])

        # Since the calculated fluxes are specific rates, we multiply them by the
        # biomass concentration to get the bulk exchange rates.
        fluxes = lex_constraints.values
        fluxes *= biomass

        # This implementation is **not** efficient, so I display the current
        # simulation time using a progress bar.
        if dynamic_system.pbar is not None:
            dynamic_system.pbar.update(1)
            dynamic_system.pbar.set_description('t = {:.3f}'.format(t))

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

        return feasibility - epsilon

    return infeasible_event


def experiment(model, exchange_rxn, y0):
    T = np.linspace(0, 1, 10)
    # Biomass, carbon source
    # y0 = [0.1, 10]

    # TODO: Extremely arbitrary!!!
    K_M = 5

    with tqdm() as pbar:
        sol = solve_ivp(
            fun=make_dynamic_system(
                model, "RPOM_provisional_biomass", exchange_rxn, K_M, pbar),
            events=[make_infeasible_event(model, exchange_rxn, K_M)],
            t_span=(T.min(), T.max()),
            y0=y0,
            t_eval=T,
            rtol=1e-3,  # 1e-6,
            atol=1e-6,  # 1e-8,
            method='RK23'
        )

    return sol


def plot_data(ivp_solution, carbon_source, ax):
    ax.plot(ivp_solution.t, ivp_solution.y.T[:, 0])
    ax2 = plt.twinx(ax)
    ax2.plot(ivp_solution.t, ivp_solution.y.T[:, 1], color='r')

    ax.set_ylabel('Biomass', color='b')
    ax2.set_ylabel(carbon_source, color='r')

    return ax, ax2


def main():
    model = read_sbml_model(MODEL)

    # TODO: Clean this up and use correct media
    min_medium = {
        "EX_ca2": 0.001692,
        "EX_cl": 0.001692,
        "EX_cobalt2": 0.001128,
        "EX_cu2": 0.001128,
        "EX_fe2": 0.002619,
        "EX_fe3": 0.002538,
        "EX_k": 0.063444,
        "EX_mg2": 0.002820,
        "EX_mn2": 0.001128,
        "EX_mobd": 0.001128,
        "EX_nh4": 3.810975,
        "EX_o2": 80.262691,
        "EX_pi": 0.722799,
        "EX_so4": 0.069699,
        "EX_thm": 0.000081,
        "EX_zn2": 0.001128
    }

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
        # Get Metabolite object and exchange reaction for the given carbon source
        metabolite = model.metabolites.query(
            lambda m: m.id == carbon_source_id)[0]
        exchange_rxn = model.exchanges.query(
            lambda r: metabolite in r.reactants)[0]
        lower_bound = float(exchange_rxn._annotation["Experimental rate"])

        model.medium = {**min_medium,
                        exchange_rxn.id: -lower_bound}

        # Run the experiment, get data
        # Biomass, carbon source
        initial = uptake_data[uptake_data["Compound"] ==
                              carbon_source]["InitialMetabolite_mM"].values[0]
        y0 = [0.1, initial]
        data = experiment(model, exchange_rxn.id, y0)

        # Plot output
        fig, ax = plt.subplots()

        # Plot depletion of carbon source over 24 hrs
        _, ax2 = plot_data(data, carbon_source_id, ax)

        # Plot fitted curve, points (?)

        rate = uptake_rates[carbon_source]
        t = np.linspace(0, 1, 10)
        ax2.plot(t, initial * np.exp(rate * t), "k--")

        # Save figure
        fig.set_size_inches(8, 6)
        fig.tight_layout()
        fig.savefig(os.path.join(OUTDIR, f"{carbon_source} dFBA.png"))


if __name__ == "__main__":
    main()
