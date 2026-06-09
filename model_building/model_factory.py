import json
import os
import pandas as pd
import numpy as np

from argparse import ArgumentParser

from cobra.io import read_sbml_model, write_sbml_model, save_json_model
from cobra.flux_analysis import pfba
from model_building.stages import STAGE_REGISTRY


# step 1: run with nobiocyc
# step 2: regenerate biocyc templates (model_building/biocyc_update_pipeline/transform/build_templates.py)
# step 3: run with biocyc 

# DEFAULT_CONFIG = "model_building/blueprints/Rpom_05_hwa__no_biocyc.json"
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
            
            # Save spreadsheet version of model
            reactions_df, metabolites_df, genes_df = model_to_dfs(model)
            with pd.ExcelWriter(out.replace(".xml", ".xlsx")) as writer:
                reactions_df.to_excel(writer, sheet_name="Reactions", index=False)
                metabolites_df.to_excel(writer, sheet_name="Metabolites", index=False)
                genes_df.to_excel(writer, sheet_name="Genes", index=False)

        return model


def model_to_dfs(model):
    with model:
        ex_glc = model.reactions.get_by_id("EX_glc")
        ex_glc.lower_bound = -5.44
        sol_glc = pfba(model)
    
    with model:
        ex_ac = model.reactions.get_by_id("EX_ac")
        ex_ac.lower_bound = -15
        sol_ac = pfba(model)

    reactions_df = pd.DataFrame([
        {
            "ID": reaction.id,
            "name": reaction.name,
            "gpr": reaction.gene_reaction_rule,
            "reaction": reaction.reaction,
            "lb": reaction.lower_bound,
            "ub": reaction.upper_bound,
            "glucose normalized flux": abs(sol_glc.fluxes[reaction.id] / sol_glc.fluxes["Rpom_hwa_biomass"]),
            "glucose sign": np.sign(sol_glc.fluxes[reaction.id]),
            "acetate normalized flux": abs(sol_ac.fluxes[reaction.id] / sol_glc.fluxes["Rpom_hwa_biomass"]),
            "acetate sign": np.sign(sol_ac.fluxes[reaction.id]),
            **{f"notes__{key}" : value for key, value in reaction.notes.items()},
            **{f"annotation__{key}" : value for key, value in reaction.annotation.items()}
        }
        for reaction in model.reactions
    ])

    metabolites_df = pd.DataFrame([
        {
            "ID": metabolite.id,
            "name": metabolite.name,
            "formula": metabolite.formula,
            "compartment": metabolite.compartment,
            "charge": metabolite.charge,
            "reactions": ", ".join(rxn.id for rxn in metabolite.reactions),
            **{f"notes__{key}": value for key, value in metabolite.notes.items()},
            **{f"annotation__{key}": value for key, value in metabolite.annotation.items()}
        }
        for metabolite in model.metabolites
    ])

    genes_df = pd.DataFrame([
        {
            "ID": gene.id,
            "name": gene.name,
            "reactions": ", ".join(rxn.id for rxn in gene.reactions),
            **{f"notes__{key}": value for key, value in gene.notes.items()},
            **{f"annotation__{key}": value for key, value in gene.annotation.items()}
        }
        for gene in model.genes
    ])

    return reactions_df, metabolites_df, genes_df


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
