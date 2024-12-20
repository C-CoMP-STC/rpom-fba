import abc
import json
import pickle
from pathlib import Path

import git
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
            raise ValueError("AddReactions stage requires a list of files/folders with reactions to add.")
        
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
            raise ValueError("RemoveReactions stage requires a list of files/folders with reactions to add.")
        
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
            raise ValueError("ModifyReactions stage requires a list of files/folders with reactions to add.")
        
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
                    model.metabolites.get_by_id(met) : coeff
                    for met, coeff in reaction["metabolites"].items()
                }
                rxn.subtract_metabolites(rxn.metabolites)
                rxn.add_metabolites(metabolites)
            
            if "annotation" in reaction:
                rxn.annotation.update(reaction["annotation"])

            rxn.gene_reaction_rule = reaction.get("gene_reaction_rule", rxn.gene_reaction_rule)

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
            raise ValueError("AnnotateReactions stage requires a path to a .pkl file, containing a dictionary of reaction annotations.")
        
        with open(params, "rb") as f:
            annotations = pickle.load(f)

        stems = annotations["stems"]
        pathways = annotations["pathways"]

        for reaction in model.reactions:
            reaction.annotation["stem"] = stems.get(reaction.id, "")
            reaction.annotation["pathways"] = list(pathways.get(reaction.id, []))
        
        return model