import abc
import json
import pickle
import warnings
import networkx as nx
import pandas as pd
from pathlib import Path

import git
from cobra.manipulation import prune_unused_metabolites, prune_unused_reactions
from cobra.manipulation.delete import remove_genes
from cobra.manipulation.validate import check_mass_balance
from cobra.core.metabolite import Metabolite
from cobra.core.model import Model
from cobra.core.reaction import Reaction
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import read_sbml_model
from macaw.main import dead_end_test, dilution_test, diphosphate_test, duplicate_test, loop_test
from pyvis.network import Network

from model_building.biomass import (add_ecoli_core_biomass_to_model,
                                    add_ecoli_full_biomass_to_model,
                                    add_hwa_biomass_to_model)
from model_building.metabolites.metabolites import ADDED_METABOLITES
from model_building.uptake_rxns import add_uptake_reactions, get_uptake_data
from utils.cobra_utils import get_or_create_exchange, set_active_bound
from utils.network import model_to_network

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
        
        model = read_sbml_model(params)
        print("Loaded base model.")

        return model


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

        # Print summary of changes
        print(f"Set medium with {len(final_medium)} exchange reactions.")

        return model


@register_stage
class AddMetabolites(Stage):
    def process(self, model: Model, params: object) -> Model:
        model.add_metabolites(ADDED_METABOLITES)
        return model


@register_stage
class RemoveMetabolites(Stage):
    def process(self, model, params):
        if not isinstance(params, list):
            raise ValueError(
                "RemoveMetabolites stage requires a list of files/folders with reactions to add.")

        root = Path(r".")

        # Collect metabolites from all files
        mets_to_remove = []
        for path in params:
            for filepath in root.glob(path):
                with open(filepath, "r") as f:
                    met_ids = json.load(f)
                    for met_id in met_ids:
                        try:
                            met = model.metabolites.get_by_id(met_id)
                            mets_to_remove.append(met)
                        except KeyError:
                            warnings.warn(f"Cannot remove {met_id}, since it is not present in the model.")

        # Remove metabolites AND associated reactions (destructive=True).
        n_reactions = len(model.reactions)
        model.remove_metabolites(mets_to_remove, destructive=True)

        print(f"Removed {len(mets_to_remove)} metabolites, and {n_reactions - len(model.reactions)} associated reactions.")

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

        # reactions = []
        for reaction in reactions_to_add:
            # Allows the use of strings as comments
            if not isinstance(reaction, dict):
                continue

            # Create a new reaction object
            rxn = Reaction()
            rxn.id = reaction["id"]
            rxn.name = reaction["name"]
            rxn.subsystem = reaction["subsystem"]

            # Skip reactions that are already in the model (they are ignored anyway by model.add_reactions,
            # but we do it here to skip subsequent steps)
            if rxn.id in model.reactions:
                print(f"Ignoring reaction '{rxn.id}' since it already exists.")
                continue

            # Add to model (need to do it now to allow building from reaction string)
            model.add_reactions([rxn])

            # reaction["metabolites"] may be a dict of metabolite ids : coefficients,
            # or a reaction string like "A + B <=> C + D"
            if isinstance(reaction["metabolites"], dict):
                metabolites = {
                    model.metabolites.get_by_id(m): v
                    for m, v in reaction["metabolites"].items()
                }
                rxn.add_metabolites(metabolites)
            elif isinstance(reaction["metabolites"], str):
                # Parse the reaction string
                reaction_str = reaction["metabolites"]
                rxn.build_reaction_from_string(reaction_str)
            else:
                raise ValueError(
                    f"Invalid format for metabolites in reaction {reaction['id']}.")

            # Set bounds and gene reaction rule
            rxn.lower_bound = reaction["lower_bound"]
            rxn.upper_bound = reaction["upper_bound"]
            rxn.gene_reaction_rule = reaction["gene_reaction_rule"]

            # Annotate reaction, gene source as manual
            rxn.annotation["source"] = "manual"
            if reaction["gene_reaction_rule"] != "":
                rxn.annotation["Gene Source"] = "manual"

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
        model.remove_reactions(reactions_to_remove, remove_orphans = True)

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

            if reaction["id"] not in model.reactions:
                warnings.warn(f"Reaction {reaction['id']} not found in model. Skipping modification.")
                continue

            rxn = model.reactions.get_by_id(reaction["id"])
            rxn.name = reaction.get("name", rxn.name)

            rxn.subsystem = reaction.get("subsystem", rxn.subsystem)
            rxn.lower_bound = reaction.get("lower_bound", rxn.lower_bound)
            rxn.upper_bound = reaction.get("upper_bound", rxn.upper_bound)

            if "metabolites" in reaction:
                # Remove current metabolites to overwrite
                rxn.subtract_metabolites(rxn.metabolites)

                # reaction["metabolites"] may be a dict of metabolite ids : coefficients,
                # or a reaction string like "A + B <=> C + D"
                if isinstance(reaction["metabolites"], dict):
                    metabolites = {
                        model.metabolites.get_by_id(m): v
                        for m, v in reaction["metabolites"].items()
                    }
                    rxn.add_metabolites(metabolites)
                elif isinstance(reaction["metabolites"], str):
                    # Parse the reaction string
                    reaction_str = reaction["metabolites"]
                    rxn.build_reaction_from_string(reaction_str)
                else:
                    raise ValueError(
                        f"Invalid format for metabolites in reaction {reaction['id']}.")
                

            if "annotation" in reaction:
                rxn.annotation.update(reaction["annotation"])

            # Update gene reaction rule if provided,
            # recording the source as manual to prevent downstream overrides
            if "gene_reaction_rule" in reaction:
                rxn.gene_reaction_rule = reaction["gene_reaction_rule"]
                rxn.annotation["Gene Source"] = "manual"

        return model


# TODO: Remove? Check if this is still necessary
@register_stage
class AddUptakeReactions(Stage):
    def process(self, model: Model, params: object) -> Model:
        # Get uptake rates and genes
        uptake_data = get_uptake_data()
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
        if not isinstance(params, dict):
            raise ValueError(
                "BioCycUpdates stage requires a dictionary of parameters.")
        
        try:
            template = params["template"]
        except KeyError:
            raise ValueError(
                "BioCycUpdates stage requires a 'template' key in the parameters dictionary, pointing to an Excel file with updates.")

        do_warnings = params.get("warnings", True)

        # Load templates for existing reactions, new reactions, new metabolites, and new genes
        existing_reactions = pd.read_excel(
            template, sheet_name="Existing reactions")
        new_reactions = pd.read_excel(template, sheet_name="New reactions")
        new_metabolites = pd.read_excel(template, sheet_name="New metabolites")
        new_genes = pd.read_excel(template, sheet_name="New genes")

        # Keep track of initial number of reactions, metabolites, and genes for later logging
        initial_reactions = len(model.reactions)
        initial_metabolites = len(model.metabolites)
        initial_genes = len(model.genes)

        # First, clear existing genes and gene rules (which use bad ids)
        # (except for manually annotated genes)
        for reaction in model.reactions:
            if reaction.annotation.get("Gene Source") != "manual":
                reaction.gene_reaction_rule = ""
        genes_to_remove = []
        for gene in model.genes:
            if len(gene.reactions) ==0:
                genes_to_remove.append(gene)
        remove_genes(model, genes_to_remove, remove_reactions=False)

        # Update existing reactions, by adding gene rules and removing unused reactions
        genes_added = set()
        removed_reactions = []
        for _, row in existing_reactions.iterrows():
            reaction = model.reactions.get_by_id(row["Reaction ID"])

            # # Get gene rule, preferring DB2 if available
            # in_db1 = row["In DB1?"]
            # in_db2 = row["In DB2?"]
            # if in_db2:
            #     gene_reaction_rule = row["Gene-reaction rule 2"]
            # elif in_db1:
            #     gene_reaction_rule = row["Gene-reaction rule 1"]
            # else:
            #     gene_reaction_rule = None

            # Get gene rule, preferring DB1 if available
            in_db1 = row["In DB1?"]
            in_db2 = row["In DB2?"]
            if in_db1:
                gene_reaction_rule = row["Gene-reaction rule 1"]
            elif in_db2:
                gene_reaction_rule = row["Gene-reaction rule 2"]
            else:
                gene_reaction_rule = None

            # Update gene rule if we have one (do not overwrite manual annotations)
            # (Also skipping empty gene rules of the form "()" since the parser doesn't like them)
            if not pd.isna(gene_reaction_rule) and gene_reaction_rule != "()":
                if reaction.annotation.get("Gene Source") != "manual":
                    reaction.gene_reaction_rule = gene_reaction_rule
                    genes_added.update(reaction.genes)
                    for gene in reaction.genes:
                        gene.annotation["source"] = ("Ruegeria pomeroyi DSS-3"
                                                    if in_db1
                                                    else "Ruegeria pomeroyi DSS-3 representative genome")

            # Remove reaction if unused and not manually added
            if (row["Recommendation"] == "Delete"
                and reaction.annotation.get("source") != "manual"):
                model.remove_reactions([reaction])
                removed_reactions.append(row["Reaction ID"])
        
        # Add new reactions
        added_reactions = []
        added_metabolites = []
        for _, row in new_reactions.iterrows():
            if row["Recommendation"] != "Add":
                continue

            # Get database of origin (preferring DB1)
            in_db1 = row["In DB1?"]
            in_db2 = row["In DB2?"]
            
            # Get bounds, defaulting to (-1000, 1000) if not found
            bounds = row["Bounds 1"] if in_db1 else row["Bounds 2"]
            if not pd.isna(bounds):
                bounds = eval(bounds)
            else:
                bounds = (-1000, 1000)
                if do_warnings:
                    warning = f"Bounds not found for reaction {row['Reaction ID']}, assuming reversible."
                    warnings.warn(warning)
            
            # Get reaction common name, if any
            name = row["Reaction name 1" if in_db1 else "Reaction name 2"]
            name = name if not pd.isna(name) else ""

            # Build reaction object
            reaction = Reaction(id=row["Reaction ID"],
                                name=name,
                                lower_bound=bounds[0],
                                upper_bound=bounds[1])
            
            # Add gene-reaction-rule (do not overwrite manual annotations)
            gene_reaction_rule = row["Gene-reaction rule 1" if in_db1 else "Gene-reaction rule 2"]
            if reaction.annotation.get("Gene Source") != "manual":
                if not pd.isna(gene_reaction_rule) and gene_reaction_rule != "()":
                    reaction.gene_reaction_rule = gene_reaction_rule
            
            # Build stoichiometry, incorporating new metabolites as needed
            stoichiometry = eval(row["Stoichiometry 1"] if in_db1 else row["Stoichiometry 2"])
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
            reaction.annotation["source"] = ("Ruegeria pomeroyi DSS-3"
                                            if in_db1
                                            else "Ruegeria pomeroyi DSS-3 representative genome")
            
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
            gene.annotation["source"] = ("Ruegeria pomeroyi DSS-3"
                                        if new_genes_row["In DB1?"]
                                        else "Ruegeria pomeroyi DSS-3 representative genome")
            
            # Position
            gene.annotation["left-end-position"] = str(new_genes_row["Left-Position"]) if not pd.isna(new_genes_row["Left-Position"]) else None
            gene.annotation["right-end-position"] = str(new_genes_row["Right-Position"]) if not pd.isna(new_genes_row["Right-Position"]) else None

            # Replicon
            gene.annotation["replicon"] = new_genes_row["Replicon"] if not pd.isna(new_genes_row["Replicon"]) else None
        
        print(f"Before updates, had {initial_reactions} reactions, {initial_metabolites} metabolites, and {initial_genes} genes.")
        print(f"- Removed {len(removed_reactions)} reactions.")
        print(f"- Removed {len(genes_to_remove)} genes.")
        print(f"+ Added {len(added_reactions)} reactions, {len(added_metabolites)} metabolites, and {len(genes_added)} genes.")
        print(f"After updates, have {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, and {len(model.genes)} genes.")

        return model


@register_stage
class DeadEnds(Stage):
    DEFAULT = {
        "outdir": "model/logs/"
    }
    def process(self, model: Model, params: object) -> Model:
        if params is None:
            params = self.DEFAULT
        elif isinstance(params, dict):
            params = {**self.DEFAULT, **params}
        else:
            raise ValueError(
                "DeadEnds stage requires a dictionary with keys 'outdir' (str) to store the logs, and 'remove' (bool) indicating whether to remove dead-ends.")

        # Dead ends ===============================================================================
        print("Dead-end test... (Ctrl+C to skip)")
        test_finished = False
        try:
            with model:
                for exchange in model.exchanges:
                    exchange.bounds = (-100, 1000)
                reaction_df, edgelist = dead_end_test(model, verbose=0)
                test_finished = True
        except KeyboardInterrupt:
            print("Dead-end test interrupted. Skipping.")

        if test_finished:
            reaction_df["pathways"] = [model.reactions.get_by_id(rxn).annotation.get("pathways", []) for rxn in reaction_df["reaction_id"]]
            reaction_df.to_csv(
                Path(params["outdir"]) / "dead_end_test.csv",
                index=False)

            g = nx.from_edgelist(edgelist)
            nt = Network()
            nt.from_nx(g)
            nt.show_buttons()
            nt.show(str(Path(params["outdir"]) / "dead_end_test.html"), notebook=False)
    
        return model
        

@register_stage
class LogMACAW(Stage):
    DEFAULT = {
        "outdir": "model/logs/",
        "substrates": ["EX_glc", "EX_ac"]
    }

    def process(self, model: Model, params: object) -> Model:
        if params is None:
            params = self.DEFAULT
        elif isinstance(params, dict):
            params = {**self.DEFAULT, **params}            
        else:
            raise ValueError(
                "LogMACAW stage requires a dictionary with keys 'outdir' to store the logs, and 'substrates' to specify the substrate exchanges.")

        # Create the directory if it doesn't exist
        Path(params["outdir"]).mkdir(parents=True, exist_ok=True)

        # Get flux limits on glucose to check if reactions are used
        with model:
            ex_glc = model.reactions.get_by_id("EX_glc")
            ex_glc.lower_bound = -10
            glc_fva = flux_variability_analysis(model)

        # Check mass balance
        print("Checking mass balance...")
        ignore = [
            "PHB-SYNTHESIS",
            "PHB-DEGRADATION",
            "Rpom_hwa_biomass",
            "DNA-RXN",
            "RNA-RXN",
            "PROTEIN-RXN",
            "LIPID-RXN",
            "MUREIN-RXN",
            "COFACTOR-RXN",
            "IONS-RXN",
            "PHB-STORAGE-RXN"
        ]
        mass_balance = {
            rxn: elems
            for rxn, elems in check_mass_balance(model).items()
            if (rxn not in model.boundary
                and rxn.id not in ignore)}
        mass_balance_df = []
        for rxn, elems in mass_balance.items():
            has_unknown_formulas = any([met.formula is None or met.formula == "" for met in rxn.metabolites])
            mass_balance_df.append({
                "ID": rxn.id,
                "reaction": rxn.reaction,
                "mass_balance": elems,
                "has_unknown_formulas": has_unknown_formulas,
                "glucose flux min": glc_fva.loc[rxn.id]["minimum"],
                "glucose flux max": glc_fva.loc[rxn.id]["maximum"]})
        mass_balance_df = pd.DataFrame(mass_balance_df)
        mass_balance_df.to_csv(
            Path(params["outdir"]) / "mass_balance.csv",
            index=False)
        print(f"Mass balance: {len(mass_balance_df)} reactions with mass balance issues.")

        # Run MACAW tests
        print("Running MACAW tests...")
        
        # Loops ===================================================================================
        print("Running loop test... (Ctrl+C to skip)")
        test_finished = False
        try:
            reaction_df, edgelist = loop_test(model, verbose=1)
            test_finished = True
        except KeyboardInterrupt:
            # Loop test takes a long time - building in a lazy escape hatch
            print("Loop test interrupted. Skipping.")

        if test_finished:
            reaction_df["pathways"] = [model.reactions.get_by_id(rxn).annotation.get("pathways", []) for rxn in reaction_df["reaction_id"]]
            reaction_df.to_csv(
                Path(params["outdir"]) / "loop_test.csv",
                index=False)
            
            g = nx.from_edgelist(edgelist)
            nt = Network()
            nt.from_nx(g)
            nt.show_buttons()
            nt.show(str(Path(params["outdir"]) / "loop_test.html"), notebook=False)

        # Diphosphate test ========================================================================
        print("Running diphosphate test... (Ctrl+C to skip)")
        test_finished = False
        try:
            reaction_df = diphosphate_test(model,
                                           ppi_ids=["PPI[c]", "PPI[p]"],
                                           pi_ids=["Pi[c]", "Pi[p]", "Pi[e]"],
                                           verbose=1)
            test_finished = True
        except KeyboardInterrupt:
            # Diphosphate test takes a long time - building in a lazy escape hatch
            print("Diphosphate test interrupted. Skipping.")
        
        if test_finished:
            reaction_df.to_csv(
                Path(params["outdir"]) / "diphosphate_test.csv",
                index=False
            )
        
        return model


@register_stage
class SanityChecks(Stage):
    def process(self, model: Model, params: object) -> Model:
        # Check growth rates on glucose, acetate
        print("Running sanity checks...")
        for substrate in ["EX_glc", "EX_ac"]:
            with model:
                ex_substrate = model.reactions.get_by_id(substrate)
                ex_substrate.lower_bound = -10
                growth_rate = model.optimize().objective_value
                
                message = f"Growth rate on 10 mmol/gDCW/hr {substrate}: {growth_rate:.4f} h^-1"
                print(f"\033[91m{message}\033[0m" if growth_rate < 1e-2 else message)
        
        # Check that biomass, ATP cannot be generated without carbon
        with model:
            # Set all carbon sources to 0
            for ex in model.exchanges:
                if "C" in list(ex.metabolites)[0].formula:
                    ex.lower_bound = 0
                    ex.upper_bound = 0
            
            # Set maintenance flux to 0
            atpm = model.reactions.get_by_id("ATPM")
            atpm.bounds = (0, 0)
            
            # Check biomass production
            with model:
                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mBiomass production is possible without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mBiomass production is not possible without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")

            # Check ATP production
            with model:
                atp_sink = model.add_boundary(
                    model.metabolites.get_by_id("ATP[c]"),
                    type="sink",
                    reaction_id="ATP-sink",
                    lb= 0,
                    ub= 1000,)
                model.objective = atp_sink
                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mATP production is possible without carbon sources!  (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mATP production is not possible without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")
        
            # Check ATPM
            with model:
                atpm.bounds = (0, 1000)
                model.objective = atpm
                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mATPM can sustain flux without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                    for rxn in sol.fluxes[sol.fluxes != 0].index:
                        print(f"\t{rxn}: {sol.fluxes[rxn]:.2f}\n\t\t{model.reactions.get_by_id(rxn).reaction}")
                else:
                    print(f"\033[92mATPM cannot sustain flux without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")
            
            # Check NADH sink
            with model:
                nadh_sink = model.add_boundary(
                    model.metabolites.get_by_id("NADH[c]"),
                    type="sink",
                    reaction_id="NADH-sink",
                    lb= 0,
                    ub= 1000,)
                model.objective = nadh_sink
                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mNADH can be produced without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mNADH cannot be produced without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")

            # Check NADPH sink
            with model:
                nadph_sink = model.add_boundary(
                    model.metabolites.get_by_id("NADPH[c]"),
                    type="sink",
                    reaction_id="NADPH-sink",
                    lb= 0,
                    ub= 1000,)
                model.objective = nadph_sink
                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mNADPH can be produced without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mNADPH cannot be produced without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")
            
            # Check FADH2 sink
            with model:
                fadh2_sink = model.add_boundary(
                    model.metabolites.get_by_id("FADH2[c]"),
                    type="sink",
                    reaction_id="FADH2-sink",
                    lb= 0,
                    ub= 1000,)
                model.objective = fadh2_sink
                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mFADH2 can be produced without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mFADH2 cannot be produced without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")


            # Check oxidizing NADH
            with model:
                nadhm = Reaction("NADHM", lower_bound=0, upper_bound=1000)
                nadhm.add_metabolites(
                    {
                        model.metabolites.get_by_id("NADH[c]"): -2,
                        model.metabolites.get_by_id("PROTON[c]") : -2,
                        model.metabolites.get_by_id("OXYGEN-MOLECULE[c]"): -1,
                        model.metabolites.get_by_id("NAD[c]"): 2,
                        model.metabolites.get_by_id("WATER[c]"): 2
                    }
                )
                model.add_reactions([nadhm])
                model.objective = nadhm

                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mNADH can be oxidized without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mNADH cannot be oxidized without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")
            
            # Check oxidizing NADPH
            with model:
                nadphm = Reaction("NADHM", lower_bound=0, upper_bound=1000)
                nadphm.add_metabolites(
                    {
                        model.metabolites.get_by_id("NADPH[c]"): -2,
                        model.metabolites.get_by_id("PROTON[c]") : -2,
                        model.metabolites.get_by_id("OXYGEN-MOLECULE[c]"): -1,
                        model.metabolites.get_by_id("NADP[c]"): 2,
                        model.metabolites.get_by_id("WATER[c]"): 2
                    }
                )
                model.add_reactions([nadphm])
                model.objective = nadphm

                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mNADPH can be oxidized without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mNADPH cannot be oxidized without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")
            
            # FADH2
            with model:
                fadh2m = Reaction("FADH2M", lower_bound = 0, upper_bound= 1000)
                fadh2m.add_metabolites(
                    {
                        model.metabolites.get_by_id("FADH2[c]"): -2,
                        model.metabolites.get_by_id("OXYGEN-MOLECULE[c]"): -1,
                        model.metabolites.get_by_id("FAD[c]"): 2,
                        model.metabolites.get_by_id("PROTON[c]"): 2,
                        model.metabolites.get_by_id("WATER[c]"): 2
                    }
                )
                model.add_reactions([fadh2m])
                model.objective = fadh2m

                sol = model.optimize()
                if sol.objective_value > 0:
                    print(f"\033[91mFADH2 can be oxidized without carbon sources! (flux = {sol.objective_value:.2f})\033[0m")
                else:
                    print(f"\033[92mFADH2 cannot be oxidized without carbon sources. (flux = {sol.objective_value:.2f})\033[0m")


        return model


@register_stage
class NetworkChecks(Stage):
    DEFAULT = {
        "outdir": "model/logs/"
    }

    def process(self, model, params) -> Model:
        """
        Conduct network-based checks (just checking number of connected components for now.)
        """
        if params is None:
            params = self.DEFAULT
        elif isinstance(params, dict):
            params = {**self.DEFAULT, **params}            
        else:
            raise ValueError(
                "NetworkChecks stage requires a dictionary with key 'outdir' to store the logs.")

        # Create the directory if it doesn't exist
        Path(params["outdir"]).mkdir(parents=True, exist_ok=True)

        # Create graph representation of model
        g = model_to_network(model)

        # Get number of connected components and warn if there are >1.
        components = list(nx.connected_components(g))

        if len(components) == 1:
            print("\033[92mNetwork has one connected component.\033[0m")
        else:
            print(f"\033[91mNetwork has {len(components)} connected components!\033[0m")
            print("\033[91mCheck logs for component memberships.\033[0m")
        
        components_log = []
        for i, component in enumerate(components):
            for elem in component:
                components_log.append(
                    {
                        "ID" : elem,
                        "Component" : i
                    }
                )
        components_log = pd.DataFrame(components_log)

        with open(Path(params["outdir"]) / "connected_components.csv", "w") as f:
            components_log.to_csv(f, index=False)

        return model


@register_stage
class RemoveOrphans(Stage):
    def process(self, model: Model, params: object) -> Model:
        """
        Remove orphan metabolites, genes, and reactions from the model.
        """
        print(f"Before removing orphans: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes.")

        model, removed_mets = prune_unused_metabolites(model)
        model, removed_reactions = prune_unused_reactions(model)
        
        orphan_genes = [gene for gene in model.genes if len(gene.reactions) == 0]
        remove_genes(model, orphan_genes)

        print(f"Removed {len(removed_mets)} orphan metabolites, {len(removed_reactions)} orphan reactions, and {len(orphan_genes)} orphan genes.")
        print(f"After removing orphans: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes.")

        return model
