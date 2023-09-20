import json
import os
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
from cobra.io import read_sbml_model

from experiments.fast_dFBA import setup_drawdown
from utils.cobra_utils import get_or_create_exchange, set_active_bound

MODEL = "model/Rpom_05.xml"
OUTDIR = "out/feasibility/o2/"
CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"
GROWTH_RATES = "parameters/growth_rates/fitted_growth_rates.csv"

def main():
    model = read_sbml_model(MODEL)
    setup_drawdown(model)

    # Load carbon sources to test
    with open(CARBON_SOURCES, "r") as f:
        carbon_sources = json.load(f)
    carbon_sources["acetate"] = "ACET[e]"

    # Load experimental growth rates
    growth_rates = pd.read_csv(GROWTH_RATES)

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    # Get relevant reactions
    ex_o2 = model.reactions.get_by_id("EX_o2")
    atpm = model.reactions.get_by_id("ATPM")

    for carbon_source, carbon_source_id in carbon_sources.items():
        with model:
            # Get exchange reaction and calculated rate for the given carbon source
            exchange_rxn = get_or_create_exchange(model, carbon_source_id)
            experimental_rate = abs(float(exchange_rxn.annotation["Experimental rate"]))
            set_active_bound(exchange_rxn, experimental_rate)

            # Sweep oxygen bound, keeping track of objective value
            # with/without maintenance flux
            x = np.linspace(0, 100, 100)
            y_atpm = np.zeros_like(x)
            y_no_atpm = np.zeros_like(x)

            # With maintenance
            for i, v in enumerate(x):
                set_active_bound(ex_o2, v)
                solution = model.optimize()
                if solution.status != "infeasible":
                    y_atpm[i] = solution.objective_value

            # Without maintenance
            atpm.bounds = (0, 0)
            for i, v in enumerate(x):
                set_active_bound(ex_o2, v)
                solution = model.optimize()
                if solution.status != "infeasible":
                    y_no_atpm[i] = solution.objective_value
                    
            fig, ax = plt.subplots()
            ax.plot(x, y_atpm, "b-", label="Maintenance")
            ax.plot(x, y_no_atpm, "b--", label="No maintenance")

            height = y_atpm.max()
            ax.fill_between(x, height * (y_atpm > 0), color="g", step="pre", alpha=0.4, label="Feasibility")
            ax.set_ylim(0, height)

            if carbon_source in growth_rates.columns:
                ax.hlines([growth_rates[carbon_source].mean()], 0, 100, ["r"], "solid", label="Experimental")
        
            # Labels
            ax.set_xlabel(f"O2 supply ($\\frac{{mmol}}{{gDCW \\cdot hr}}$)")
            ax.set_ylabel(r"Growth rate ($hr^{-1}$)")
            ax.set_title(f"O2 limitation on {carbon_source}")
            
            # Save figure
            fig.set_size_inches(4, 2)
            fig.tight_layout()
            fig.savefig(os.path.join(OUTDIR, f"{carbon_source} O2 feasibility scan.png"))


if __name__ == "__main__":
    main()
