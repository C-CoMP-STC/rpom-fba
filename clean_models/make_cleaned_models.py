"""
clean_models.py

Script for validating and cleaning the initial models (second step of the pipeline,
following model conversion from .mat to .sbml).
"""
from argparse import ArgumentParser
import os
import json
from cobra import Reaction
from cobra.io import read_sbml_model, write_sbml_model

from biomass import get_ecoli_biomass, ecoli_biomass_to_rpom, get_core_biomass_stoich
from uptake_rxns import get_uptake_rates, add_uptake_reactions


BASE_MODEL_DIR = "base_model/"
BASE_MODELS = ["Rpom_0.xml", "Rpom_02.xml", "Rpom_03.xml", "Rpom_04.xml", "Rpom_05.xml",
               "Rpom_06.xml", "Rpom_025.xml", "Rpom_035.xml", "Rpom_045.xml", "Rpom_055.xml"]
CLEAN_MODEL_DIR = "clean_models/"
MANUAL_METABOLITE_MATCHES = "clean_models/manual_met_links.json"


def main(biomass_objective):
    # Get ecoli biomass reaction stoichiometry
    ecoli_biomass = get_ecoli_biomass()

    # Get uptake rates
    uptake_rates = get_uptake_rates()

    # Load models
    # (loading all models before cleaning because of excessive warning messages)
    rpom_models = []
    for model_file in BASE_MODELS:
        model = read_sbml_model(os.path.join(BASE_MODEL_DIR, model_file))
        rpom_models.append(model)

    # Load manual metabolite id links between E. coli and R. pom.
    with open(MANUAL_METABOLITE_MATCHES, "r") as f:
        manual_matches = json.load(f)

    # Add biomass reaction, using ecoli biomass for now
    for model in rpom_models:
        match biomass_objective:
            case "ecoli":
                # Get stoichiometry of biomass reaction
                biomass_stoich = ecoli_biomass_to_rpom(
                    model, ecoli_biomass, manual_matches=manual_matches)

                # Create biomass reaction
                biomass_rxn = Reaction("RPOM_provisional_biomass",
                                    "RPOM_provisional_biomass",
                                    lower_bound=0,
                                    upper_bound=1000)
                biomass_rxn.add_metabolites(biomass_stoich)
            case "core":
                biomass_stoich = get_core_biomass_stoich()
                biomass_stoich = {model.metabolites.get_by_id(met) : coeff for met, coeff in biomass_stoich.items()}
                biomass_rxn = Reaction("RPOM_provisional_biomass",
                                       "RPOM_provisional_biomass",
                                       lower_bound=0,
                                       upper_bound=1000)
                biomass_rxn.add_metabolites(biomass_stoich)
            case _:
                raise ValueError(f"{biomass_objective} is not a recognized biomass objective.")

        # Add to model, set as objective
        model.add_reactions([biomass_rxn])
        model.objective = biomass_rxn

        # Add uptake reactions
        add_uptake_reactions(model, uptake_rates)

        # Save cleaned model
        write_sbml_model(model, os.path.join(
            CLEAN_MODEL_DIR, f"{model.id}.xml"))


if __name__ == "__main__":
    argparser = ArgumentParser("Create cleaned models from the base models.")

    argparser.add_argument("--biomass", "-b", default="core", choices=["ecoli", "core"])
    
    args = argparser.parse_args()

    main(args.biomass)
