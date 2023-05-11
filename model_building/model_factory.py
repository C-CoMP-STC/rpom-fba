import json
import os
from argparse import ArgumentParser

import git
from biomass import (add_ecoli_core_biomass_to_model,
                     add_ecoli_full_biomass_to_model)
from cobra.io import read_sbml_model, write_sbml_model
from uptake_rxns import add_uptake_reactions, get_uptake_rates


def get_git_hash():
    try:
        repo = git.Repo(search_parent_directories=True)
        sha = repo.head.object.hexsha
        return sha
    except:
        return None


class ModelFactory:
    def __init__(self, config):
        match config:
            case str():
                self.config_file = config
                with open(config, "r") as f:
                    self.config = json.load(f)
            case dict():
                self.config = config
            case _:
                raise TypeError(
                    "Invalid type for config (must be str path, or dict)")

    def build_model(self, out=None):
        # Start from chosen base model
        model = read_sbml_model(self.config["base"])

        match self.config["biomass"]:
            case "ecoli-core":
                add_ecoli_core_biomass_to_model(model)
            case "ecoli-full":
                add_ecoli_full_biomass_to_model(model)

        # Get uptake rates
        uptake_rates = get_uptake_rates()

        # Add uptake reactions
        add_uptake_reactions(model, uptake_rates)

        # Save cleaned model
        if out is not None:
            os.makedirs(os.path.dirname(out), exist_ok=True)
            write_sbml_model(model, out)

        return model


def main(config_file, out_file):
    model_factory = ModelFactory(config_file)
    model_factory.build_model(out_file)


if __name__ == "__main__":
    argparser = ArgumentParser("Create cleaned models from the base models.")

    argparser.add_argument(
        "config", 
        nargs = "?",
        default="model_building/blueprints/Rpom_05_ecoli_core.json")
    
    argparser.add_argument(
        "--out",
        "-o",
        default="model/Rpom_05.xml")

    args = argparser.parse_args()

    main(args.config, args.out)
