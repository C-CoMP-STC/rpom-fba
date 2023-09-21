import abc
import json

import git
from cobra.io import read_sbml_model
from cobra.core.model import Model
from cobra.core.metabolite import Metabolite
from cobra.core.reaction import Reaction

from biomass import (add_ecoli_core_biomass_to_model,
                     add_ecoli_full_biomass_to_model)
from uptake_rxns import add_uptake_reactions, get_uptake_data
from utils.cobra_utils import get_or_create_exchange, set_active_bound
from model_building.metabolites.metabolites import ADDED_METABOLITES

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

        return model


@register_stage
class MaintenanceFlux(Stage):
    def process(self, model: Model, params: object) -> Model:
        if not params:
            atpm = model.reactions.get_by_id("ATPM")
            atpm.bounds = (0, 0)
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
            final_medium[exchange.id] = value

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
        with open(params, "r") as f:
            reactions_to_add = json.load(f)

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
class AddUptakeReactions(Stage):
    def process(self, model: Model, params: object) -> Model:
        # Get uptake rates and genes
        uptake_data = get_uptake_data()

        # Add uptake reactions
        add_uptake_reactions(model, uptake_data)

        return model
