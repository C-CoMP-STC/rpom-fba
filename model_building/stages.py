import abc

import git
from cobra.io import read_sbml_model
from cobra.core.model import Model

from biomass import (add_ecoli_core_biomass_to_model,
                     add_ecoli_full_biomass_to_model)
from uptake_rxns import add_uptake_reactions, get_uptake_rates

STAGE_REGISTRY = {}

def register_stage(cls):
    STAGE_REGISTRY[cls.__name__] = cls
    return cls


class Stage(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def process(self, model : Model, params : object) -> Model:
        pass


@register_stage
class BaseModel(Stage):
    def process(self, model: Model, params: object) -> Model:
        if not isinstance(params, str):
            raise TypeError(f"BaseModel stage requires a path to a base model (as a str). Got {type(params)}.")

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
        # Add chosen biomass objective
        match params:
            case "ecoli-core":
                add_ecoli_core_biomass_to_model(model)
            case "ecoli-full":
                add_ecoli_full_biomass_to_model(model)
        
        return model


@register_stage
class AddUptakeReactions(Stage):
    def process(self, model: Model, params: object) -> Model:
        # Get uptake rates
        uptake_rates = get_uptake_rates()

        # Add uptake reactions
        add_uptake_reactions(model, uptake_rates)

        return model

