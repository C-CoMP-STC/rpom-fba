import json
import os

import matplotlib.pyplot as plt
import numpy as np
from cobra.io import read_sbml_model

from utils.cobra_utils import get_or_create_exchange, set_active_bound

MODEL = "model/Rpom_05.xml"
OUTDIR = "out/feasibility/"
CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"

def main():
    model = read_sbml_model(MODEL)

    supp_medium = {k : 1000. for k in model.medium.keys()}
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

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    for carbon_source, carbon_source_id in carbon_sources.items():
        with model:
            # Get exchange reaction and calculated rate for the given carbon source
            exchange_rxn = get_or_create_exchange(model, carbon_source_id)
            experimental_rate = abs(float(exchange_rxn.annotation["Experimental rate"]))

            # Sweep bound on exchange reaction,
            # keeping track of objective value and feasibility
            x = np.linspace(0, max(100, experimental_rate), 100)
            y = np.zeros_like(x)
            feasible = np.zeros_like(x)
            for i, v in enumerate(x):
                set_active_bound(exchange_rxn, v)
                solution = model.optimize()
                
                if solution.status != "infeasible":
                    feasible[i] = 1
                    y[i] = solution.objective_value
                    
            fig, ax = plt.subplots()
            ax.plot(x, y, "b-", label="Objective value")

            height = y.max()
            ax.fill_between(x, height * feasible, color="g", step="pre", alpha=0.4, label="Feasibility")
            ax.set_ylim(0, height)

            ax.vlines([experimental_rate], 0, max(height, 0.01), ["r"])
            
            # Labels
            ax.set_xlabel(f"{carbon_source} supply ($\\frac{{mmol}}{{gDCW \\cdot hr}}$)")
            ax.set_ylabel("Biomass flux")
            ax.set_title(f"Feasibility on {carbon_source} limitation")
            
            # Save figure
            fig.set_size_inches(4, 2)
            fig.tight_layout()
            fig.savefig(os.path.join(OUTDIR, f"{carbon_source} feasibility scan.png"))


if __name__ == "__main__":
    main()
