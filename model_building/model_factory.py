import json
import os
from argparse import ArgumentParser

from cobra.io import read_sbml_model, write_sbml_model, save_json_model
from model_building.stages import STAGE_REGISTRY


DEFAULT_CONFIG = "model_building/blueprints/Rpom_05_hwa.json"
DEFAULT_MODEL = "model/Rpom_05.xml"


class ModelFactory:
    def __init__(self, config):
        match config:
            case str():
                self.config_file = config
                with open(config, "r") as f:
                    self.config = json.load(f)
            case dict():
                self.config_file = None
                self.config = config
            case _:
                raise TypeError(
                    "Invalid type for config (must be str path, or dict)")

    def build_model(self, out=None, verbose=True):
        model = None
        for stage, params in self.config.items():
            if verbose:
                print(f"\n\033[93mRunning stage {stage} ==================================\033[0m")
                print(f"\033[36mParams: {params}\033[0m\n")
            
            # Log model stats before stage
            n_reactions_before = len(model.reactions) if model else 0
            n_metabolites_before = len(model.metabolites) if model else 0
            n_genes_before = len(model.genes) if model else 0
            
            # Run stage
            model = STAGE_REGISTRY[stage]().process(model, params)

            # Log model stats after stage
            n_reactions_after = len(model.reactions) if model else 0
            n_metabolites_after = len(model.metabolites) if model else 0
            n_genes_after = len(model.genes) if model else 0

            # Print model stats
            if verbose:
                print(f"\033[90mReactions: {n_reactions_before} -> {n_reactions_after} ({n_reactions_after - n_reactions_before})\033[0m")
                print(f"\033[90mMetabolites: {n_metabolites_before} -> {n_metabolites_after} ({n_metabolites_after - n_metabolites_before})\033[0m")
                print(f"\033[90mGenes: {n_genes_before} -> {n_genes_after} ({n_genes_after - n_genes_before})\033[0m")

        # Save cleaned model
        if out is not None:
            os.makedirs(os.path.dirname(out), exist_ok=True)
            write_sbml_model(model, out)
            save_json_model(model, out.replace(".xml", ".json"))

        return model


def rebuild_and_get_model(config_file=DEFAULT_CONFIG, model_out=DEFAULT_MODEL):
    ModelFactory(config_file).build_model(model_out)
    return read_sbml_model(model_out)


def main(config_file, out_file):
    model_factory = ModelFactory(config_file)
    model_factory.build_model(out_file)


if __name__ == "__main__":
    argparser = ArgumentParser("Create cleaned models from the base models.")

    argparser.add_argument(
        "config",
        nargs="?",
        default=DEFAULT_CONFIG)

    argparser.add_argument(
        "--out",
        "-o",
        default=DEFAULT_MODEL)

    args = argparser.parse_args()

    main(args.config, args.out)
