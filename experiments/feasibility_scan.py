import json
import os
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from cobra.io import read_sbml_model

from utils.cobra_utils import get_or_create_exchange, set_active_bound
from model_building.model_factory import rebuild_and_get_model

OUTDIR = "out/feasibility/"
CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"
GROWTH_RATES = "parameters/growth_rates/fitted_growth_rates.csv"


def main():
    model = rebuild_and_get_model()  # read_sbml_model(DEFAULT_MODEL)

    # Load carbon sources to test
    with open(CARBON_SOURCES, "r") as f:
        carbon_sources = json.load(f)
    carbon_sources["acetate"] = "ACET[e]"

    # Load expected growth rates
    growth_rates = pd.read_csv(GROWTH_RATES)

    # Get maintenance reaction
    atpm = model.reactions.get_by_id("ATPM")
    atpm_bounds = atpm.bounds

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    for carbon_source, carbon_source_id in tqdm(carbon_sources.items()):
        with model:
            # Get exchange reaction and calculated rate for the given carbon source
            exchange_rxn = get_or_create_exchange(model, carbon_source_id)
            experimental_rate = abs(
                float(exchange_rxn.annotation["Experimental rate"]))

            growth_rate = (growth_rates[carbon_source].mean()
                           if carbon_source in growth_rates.columns
                           else 0)

            # Sweep bound on exchange reaction,
            # keeping track of objective value and feasibility
            x = np.linspace(0, max(50, experimental_rate), 50)
            y_atpm = np.zeros_like(x)
            y_no_atpm = np.zeros_like(x)
            feasible = np.zeros_like(x)
            for i, v in enumerate(x):
                set_active_bound(exchange_rxn, v)
                solution = model.optimize()

                if solution.status != "infeasible":
                    feasible[i] = 1
                    y_atpm[i] = solution.objective_value

                # atpm.bounds = (0, 0)
                # solution = model.optimize()
                # atpm.bounds = atpm_bounds

                # if solution.status != "infeasible":
                #     y_no_atpm[i] = solution.objective_value

            fig, ax = plt.subplots()
            ax.plot(x, y_atpm, "b-", label="Objective value")
            ax.plot(x, y_no_atpm, "b--")

            height = y_atpm.max()
            ax.fill_between(x, height * feasible, color="g",
                            step="pre", alpha=0.4, label="Feasibility")
            ax.set_ylim(0, height)

            ax.hlines([growth_rate], 0, x.max(), ["r"])
            ax.vlines([experimental_rate], 0, max(height, 0.01), ["r"])

            # Labels
            ax.set_xlabel(
                f"{carbon_source} supply ($\\frac{{mmol}}{{gDCW \\cdot hr}}$)")
            ax.set_ylabel("Biomass flux")
            ax.set_title(f"Feasibility on {carbon_source} limitation")

            # Save figure
            fig.set_size_inches(4, 2)
            fig.tight_layout()
            fig.savefig(os.path.join(
                OUTDIR, f"{carbon_source} feasibility scan.png"))


if __name__ == "__main__":
    main()
