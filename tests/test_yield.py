from utils.units import u
from utils.cobra_utils import get_or_create_exchange, get_active_bound
from experiments.fast_dFBA import MichaelisMentenBounds
import cobra
from cobra.io import read_sbml_model
import matplotlib.pyplot as plt
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")


LOGISTIC_RATES_FILE = "parameters/logistic_fits/fitted_rates.json"
CARBON_SOURCE_IDS = "parameters/uptake_rates/carbon_source_ids.json"
YIELD_OUTDIR = "out/tests/yield/"
MODEL = "model/Rpom_05.xml"


def main():
    K_M = 1

    # Load and set up model
    model = read_sbml_model(MODEL)

    supp_medium = {k: 1000. for k in model.medium.keys()}
    supp_medium["EX_fe2"] = 1000.
    model.medium = supp_medium

    # Remove biotin from objective temporarily as biotin is blocking
    # TODO: fix biotin production?
    biotin = model.metabolites.get_by_id("BIOTIN[c]")
    biomass = model.reactions.get_by_id("RPOM_provisional_biomass")
    biomass.subtract_metabolites({biotin: biomass.metabolites[biotin]})

    with open(LOGISTIC_RATES_FILE, "r") as f:
        rates = json.load(f)

    with open(CARBON_SOURCE_IDS, "r") as f:
        carbon_source_ids = json.load(f)
    carbon_source_ids["acetate"] = "ACET[e]"

    # Make growth rate vs. [S] plots
    S_range = np.linspace(0, 10, 100) * u.mM
    for substrate, (k, alpha) in rates.items():
        mu_fit = k * S_range

        substrate_id = carbon_source_ids[substrate]
        ex = get_or_create_exchange(model, substrate_id)
        V_max = abs(float(ex._annotation["Experimental rate"]))

        model_yield = np.zeros_like(S_range)
        uptake_rate = np.zeros_like(S_range)
        for i, s in enumerate(S_range):
            with model:
                MichaelisMentenBounds(substrate_id, V_max, K_M).bound(ex, s)

                try:
                    lex_constraints = cobra.util.add_lexicographic_constraints(
                        model, ["RPOM_provisional_biomass", ex.id], ['max', "max"])
                    fluxes = lex_constraints.values
                    mu_hat = fluxes[0] # model.optimize()
                    model_yield[i] = mu_hat# .objective_value

                    # Check that more C is not being created in biomass than is being put in
                    carbon_flux_in = abs(ex.flux) * model.metabolites.get_by_id(substrate_id).elements.get("C", 0)
                    carbon_flux_through_biomass = sum(
                        [- met.elements.get("C", 0) * (mu_hat / coeff)
                            for met, coeff in biomass.metabolites.items()])

                    print(f"{substrate_id} : {carbon_flux_in} > {carbon_flux_through_biomass} ({carbon_flux_in > carbon_flux_through_biomass})")
                except Exception as e:
                    pass
                uptake_rate[i] = get_active_bound(ex)            

        fig, ax = plt.subplots()
        ax.plot(S_range, mu_fit.magnitude, "k", label="Logistic fit")
        ax.plot(S_range, model_yield, "r--", label="FBA")
        ax.legend()
        ax.set_xlabel("[S] ($mM$)")
        ax.set_ylabel("Growth rate ($hr^{-1}$)")
        ax.set_title(f"{substrate} yield")

        # ax2 = ax.twinx()
        # ax2.plot(S_range, uptake_rate, "g--")
        # ax2.set_ylabel("M-M uptake rate (mmol/gDCW/hr)")

        os.makedirs(YIELD_OUTDIR, exist_ok=True)
        fig.set_size_inches(4, 3)
        fig.tight_layout()
        fig.savefig(os.path.join(YIELD_OUTDIR, f"{substrate}_yield.png"))


if __name__ == "__main__":
    main()
