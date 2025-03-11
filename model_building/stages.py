import abc
import json
import pickle
import warnings
import pandas as pd
from pathlib import Path

import git
from cobra.manipulation.delete import remove_genes
from cobra.core.metabolite import Metabolite
from cobra.core.model import Model
from cobra.core.reaction import Reaction
from cobra.io import read_sbml_model

from model_building.biomass import (add_ecoli_core_biomass_to_model,
                                    add_ecoli_full_biomass_to_model,
                                    add_hwa_biomass_to_model)
from model_building.metabolites.metabolites import ADDED_METABOLITES
from model_building.uptake_rxns import add_uptake_reactions, get_uptake_data
from utils.cobra_utils import get_or_create_exchange, set_active_bound

STAGE_REGISTRY = {}


def register_stage(cls):
    STAGE_REGISTRY[cls.__name__] = cls
    return cls


class Stage(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def process(self, model: Model, params: object) -> Model:
        pass


@register_stage
class BaseModel(Stage):
    def process(self, model: Model, params: object) -> Model:
        if not isinstance(params, str):
            raise TypeError(
                f"BaseModel stage requires a path to a base model (as a str). Got {type(params)}.")
        return read_sbml_model(params)


@register_stage
class MetaData(Stage):
    def process(self, model: Model, params: object) -> Model:
        # Store config and git hash in annotation
        # model.annotation["config-file"] = self.config_file
        # model.annotation["config"] = json.dumps(self.config)
        model.annotation["git-hash"] = self.get_git_hash()

        return model

    def get_git_hash(self):
        try:
            repo = git.Repo(search_parent_directories=True)
            sha = repo.head.object.hexsha
            return sha
        except:
            return None


@register_stage
class BiomassObjective(Stage):
    def process(self, model: Model, params: object) -> Model:
        # Ignore previous biomass objective
        # TODO: compare this with what we have?
        model.reactions.get_by_id("BiomassRxn").bounds = (0, 0)

        # Add chosen biomass objective
        match params:
            case "ecoli-core":
                add_ecoli_core_biomass_to_model(model)
            case "ecoli-full":
                add_ecoli_full_biomass_to_model(model)
            case "hwa":
                add_ecoli_full_biomass_to_model(model)  # Debugging
                add_hwa_biomass_to_model(model)

        return model


@register_stage
class MaintenanceFlux(Stage):
    def process(self, model: Model, params: float) -> Model:
        flux = abs(params)

        atpm = model.reactions.get_by_id("ATPM")
        atpm.bounds = (flux, flux)
        return model


@register_stage
class SetMedium(Stage):
    def process(self, model: Model, params: object) -> Model:
        # reset existing medium
        # model.medium = {rxn: 0 for rxn in model.medium}

        # Load medium
        with open(params, "r") as f:
            medium = json.load(f)

        # For given medium, find or create exchange reaction for each metabolite
        final_medium = {}
        for metabolite, value in medium.items():
            exchange = get_or_create_exchange(model, metabolite, verbose=True)

            # final_medium[exchange.id] = value

            # Currently, the values are not accurate - let everything
            # be essentially unbounded except for oxygen
            if metabolite != "OXYGEN-MOLECULE[e]":
                final_medium[exchange.id] = 1000

        # Set oxygen to 20
        final_medium["EX_o2"] = 20

        model.medium = final_medium

        return model


@register_stage
class AddMetabolites(Stage):
    def process(self, model: Model, params: object) -> Model:
        model.add_metabolites(ADDED_METABOLITES)

        return model


@register_stage
class AddReactions(Stage):
    def process(self, model: Model, params: object) -> Model:
        if not isinstance(params, list):
            raise ValueError(
                "AddReactions stage requires a list of files/folders with reactions to add.")

        root = Path(r".")

        # Collect reactions from all files
        reactions_to_add = []
        for path in params:
            for filepath in root.glob(path):
                with open(filepath, "r") as f:
                    reactions_to_add += json.load(f)

        reactions = []
        for reaction in reactions_to_add:
            # Allows the use of strings as comments
            if not isinstance(reaction, dict):
                continue

            reaction["metabolites"] = {
                model.metabolites.get_by_id(m): v
                for m, v in reaction["metabolites"].items()
            }

            rxn = Reaction()
            rxn.id = reaction["id"]
            rxn.name = reaction["name"]
            rxn.subsystem = reaction["subsystem"]
            rxn.lower_bound = reaction["lower_bound"]
            rxn.upper_bound = reaction["upper_bound"]
            rxn.add_metabolites(reaction["metabolites"])
            rxn.gene_reaction_rule = reaction["gene_reaction_rule"]

            reactions.append(rxn)

        model.add_reactions(reactions)

        return model


@register_stage
class RemoveReactions(Stage):
    def process(self, model: Model, params: object) -> Model:
        if not isinstance(params, list):
            raise ValueError(
                "RemoveReactions stage requires a list of files/folders with reactions to add.")

        root = Path(r".")

        # Collect reactions from all files
        reactions_to_remove = []
        for path in params:
            for filepath in root.glob(path):
                with open(filepath, "r") as f:
                    reactions_to_remove += json.load(f)

        # Remove reactions
        model.remove_reactions(reactions_to_remove)

        return model


@register_stage
class ModifyReactions(Stage):
    def process(self, model: Model, params: str) -> Model:
        if not isinstance(params, list):
            raise ValueError(
                "ModifyReactions stage requires a list of files/folders with reactions to add.")

        root = Path(r".")

        # Collect reactions to modify from all files
        reactions_to_change = []
        for path in params:
            for filepath in root.glob(path):
                with open(filepath, "r") as f:
                    reactions_to_change += json.load(f)

        for reaction in reactions_to_change:
            # Allows the use of strings as comments:
            if not isinstance(reaction, dict):
                continue

            rxn = model.reactions.get_by_id(reaction["id"])
            rxn.name = reaction.get("name", rxn.name)

            rxn.subsystem = reaction.get("subsystem", rxn.subsystem)
            rxn.lower_bound = reaction.get("lower_bound", rxn.lower_bound)
            rxn.upper_bound = reaction.get("upper_bound", rxn.upper_bound)

            if "metabolites" in reaction:
                metabolites = {
                    model.metabolites.get_by_id(met): coeff
                    for met, coeff in reaction["metabolites"].items()
                }
                rxn.subtract_metabolites(rxn.metabolites)
                rxn.add_metabolites(metabolites)

            if "annotation" in reaction:
                rxn.annotation.update(reaction["annotation"])

            rxn.gene_reaction_rule = reaction.get(
                "gene_reaction_rule", rxn.gene_reaction_rule)

        return model


@register_stage
class AddUptakeReactions(Stage):
    def process(self, model: Model, params: object) -> Model:
        # Get uptake rates and genes
        uptake_data = get_uptake_data()

        # Add uptake reactions
        add_uptake_reactions(model, uptake_data)

        return model


@register_stage
class AnnotateReactions(Stage):
    def process(self, model: Model, params: object) -> Model:
        if params is None:
            return model

        if not isinstance(params, str):
            raise ValueError(
                "AnnotateReactions stage requires a path to a .pkl file, containing a dictionary of reaction annotations.")

        with open(params, "rb") as f:
            annotations = pickle.load(f)

        stems = annotations["stems"]
        pathways = annotations["pathways"]

        for reaction in model.reactions:
            reaction.annotation["stem"] = stems.get(reaction.id, "")
            reaction.annotation["pathways"] = list(
                pathways.get(reaction.id, []))

        return model


@register_stage
class BioCycUpdates(Stage):
    def process(self, model: Model, params: object) -> Model:
        if not isinstance(params, str):
            raise ValueError(
                "BioCycUpdates stage requires a path to a template, built from our BioCyc ETL pipeline.")

        # Load templates for existing reactions, new reactions, new metabolites, and new genes
        existing_reactions = pd.read_excel(
            params, sheet_name="Existing reactions")
        new_reactions = pd.read_excel(params, sheet_name="New reactions")
        new_metabolites = pd.read_excel(params, sheet_name="New metabolites")
        new_genes = pd.read_excel(params, sheet_name="New genes")

        # Keep track of initial number of reactions, metabolites, and genes for later logging
        initial_reactions = len(model.reactions)
        initial_metabolites = len(model.metabolites)
        initial_genes = len(model.genes)

        # First, clear existing genes and gene rules (which use bad ids)
        for reaction in model.reactions:
            reaction.gene_reaction_rule = ""
        gene_ids = set(gene.id for gene in model.genes)
        remove_genes(model, gene_ids, remove_reactions=False)

        # Update existing reactions, by adding gene rules and removing unused reactions
        genes_added = set()
        removed_reactions = []
        for _, row in existing_reactions.iterrows():
            reaction = model.reactions.get_by_id(row["Reaction ID"])

            # Get gene rule, preferring DB2 if available
            in_db1 = row["In DB1?"]
            in_db2 = row["In DB2?"]
            if in_db2:
                gene_reaction_rule = row["Gene-reaction rule 2"]
            elif in_db1:
                gene_reaction_rule = row["Gene-reaction rule 1"]
            else:
                gene_reaction_rule = None

            # Update gene rule if we have one
            # (Also skipping empty gene rules of the form "()" since the parser doesn't like them)
            if not pd.isna(gene_reaction_rule) and gene_reaction_rule != "()":
                reaction.gene_reaction_rule = gene_reaction_rule
                genes_added.update(reaction.genes)
                for gene in reaction.genes:
                    gene.annotation["source"] = ("Ruegeria pomeroyi DSS-3 representative genome"
                                                if in_db2
                                                else "Ruegeria pomeroyi DSS-3")

            # Remove reaction if unused
            if row["Recommendation"] == "Delete":
                model.remove_reactions([reaction])
                removed_reactions.append(row["Reaction ID"])
        
        # Add new reactions
        added_reactions = []
        added_metabolites = []
        for _, row in new_reactions.iterrows():
            if row["Recommendation"] != "Add":
                continue

            # Get database of origin (preferring DB2)
            in_db1 = row["In DB1?"]
            in_db2 = row["In DB2?"]
            
            # Get bounds, defaulting to (-1000, 1000) if not found
            bounds = row["Bounds 2"] if in_db2 else row["Bounds 1"]
            if not pd.isna(bounds):
                bounds = eval(bounds)
            else:
                bounds = (-1000, 1000)
                warning = f"Bounds not found for reaction {row['Reaction ID']}, assuming reversible."
                warnings.warn(warning)
            
            # Get reaction common name, if any
            name = row["Reaction name 2" if in_db2 else "Reaction name 1"]
            name = name if not pd.isna(name) else ""

            # Build reaction object
            reaction = Reaction(id=row["Reaction ID"],
                                name=name,
                                lower_bound=bounds[0],
                                upper_bound=bounds[1])
            
            # Add gene-reaction-rule
            gene_reaction_rule = row["Gene-reaction rule 2" if in_db2 else "Gene-reaction rule 1"]
            if not pd.isna(gene_reaction_rule) and gene_reaction_rule != "()":
                reaction.gene_reaction_rule = gene_reaction_rule
            
            # Build stoichiometry, incorporating new metabolites as needed
            stoichiometry = eval(row["Stoichiometry 2"] if in_db2 else row["Stoichiometry 1"])
            metabolites = {}
            for met_id, coeff in stoichiometry.items():
                if met_id not in model.metabolites:
                    met_row = new_metabolites[new_metabolites["Metabolite ID"] == met_id.split('[')[0]].iloc[0]
                    met = Metabolite(id=met_id,
                                    name=met_row["Metabolite Name"] if not pd.isna(met_row["Metabolite Name"]) else None,
                                    formula=met_row["Formula"] if not pd.isna(met_row["Formula"]) else None,
                                    charge=met_row["Charge"] if not pd.isna(met_row["Charge"]) else None,
                                    compartment=met_id.split('[')[1][:-1])
                    model.add_metabolites([met])
                    added_metabolites.append(met_id)
                metabolites[model.metabolites.get_by_id(met_id)] = coeff
            
            # Add stoichiometry to reaction
            reaction.add_metabolites(metabolites)

            # Keep track of which database the reaction came from
            reaction.annotation["source"] = ("Ruegeria pomeroyi DSS-3 representative genome"
                                            if in_db2
                                            else "Ruegeria pomeroyi DSS-3")
            
            # Add reaction to model, keep track of new genes
            model.add_reactions([reaction])
            genes_added.update(reaction.genes)
            added_reactions.append(reaction.id)

        # Add Annotations on new genes
        for gene in genes_added:
            new_genes_row = new_genes[new_genes["Gene ID"] == gene.id].iloc[0]
            
            # Name, synonyms
            gene.name = new_genes_row["Gene Name"] if not pd.isna(new_genes_row["Gene Name"]) else ""
            gene.annotation["synonyms"] = new_genes_row["Synonyms"] if not pd.isna(new_genes_row["Synonyms"]) else []

            # database of origin
            gene.annotation["source"] = ("Ruegeria pomeroyi DSS-3 representative genome"
                                        if new_genes_row["In DB2?"]
                                        else "Ruegeria pomeroyi DSS-3")
            
            # Position
            gene.annotation["left-end-position"] = str(new_genes_row["Left-Position"]) if not pd.isna(new_genes_row["Left-Position"]) else None
            gene.annotation["right-end-position"] = str(new_genes_row["Right-Position"]) if not pd.isna(new_genes_row["Right-Position"]) else None

            # Replicon
            gene.annotation["replicon"] = new_genes_row["Replicon"] if not pd.isna(new_genes_row["Replicon"]) else None
        
        print(f"Before updates, had {initial_reactions} reactions, {initial_metabolites} metabolites, and {initial_genes} genes.")
        print(f"Removed {len(removed_reactions)} reactions.")
        print(f"Added {len(added_reactions)} reactions, {len(added_metabolites)} metabolites, and {len(genes_added)} genes.")

        return model
