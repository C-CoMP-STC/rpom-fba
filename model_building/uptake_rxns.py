import json
import re

from cobra.core.metabolite import Metabolite
from cobra.core.reaction import Reaction

from utils.cobra_utils import change_compartment, get_or_create_exchange

UPTAKE_RATES_DATA = "parameters/uptake_rates/fitted_uptake_rates.json"
UPTAKE_GENES = "parameters/uptake_rates/uptake_genes.json"
CARBON_SOURCE_IDS = "parameters/uptake_rates/carbon_source_ids.json"


def get_uptake_data(uptake_rates_file=UPTAKE_RATES_DATA,
                    uptake_genes_file=UPTAKE_GENES,
                    carbon_sources_file=CARBON_SOURCE_IDS):
    with open(uptake_rates_file, "r") as f:
        uptake_rates = json.load(f)

    with open(uptake_genes_file, "r") as f:
        uptake_genes = json.load(f)

    with open(carbon_sources_file, "r") as f:
        carbon_sources = json.load(f)
    carbon_sources["acetate"] = "ACET[e]"
    
    result = {
        carbon_sources[metabolite]: {
            "rate": rate,
            "gene_reaction_rule": uptake_genes[metabolite]
        }
        for metabolite, rate in uptake_rates.items()
    }

    return result


def add_uptake_reactions(model, uptake_data):
    # Add uptake reactions one-by-one
    for met_id, data in uptake_data.items():
        rate = data["rate"]
        gene_rule = data["gene_reaction_rule"]

        # Check if exchange reaction already exists,
        # create exchange reaction if it does not
        try:
            exchange = get_or_create_exchange(model, met_id, verbose=True)
        except KeyError:
            print(f"{met_id} not found in model, failed to create exchange reaction.")
            continue

        # Check if transport reaction already exists,
        # create transport reaction if it does not
        met_ex = model.metabolites.get_by_id(met_id)
        met_c = model.metabolites.get_by_id(change_compartment(met_id, "c"))

        transport_rxns = {
            rxn for rxn in model.reactions
            if met_ex in rxn.metabolites and met_c in rxn.metabolites
        }

        # Add reactions to transport to periplasm, if they exist
        try:
            met_p = model.metabolites.get_by_id(change_compartment(met_id, "p"))
            transport_rxns.update(rxn for rxn in model.reactions
                                  if met_ex in rxn.metabolites and met_p in rxn.metabolites)
        except:
            pass

        if len(transport_rxns) > 0:
            print(f"Found {len(transport_rxns)} transport reactions for {met_id} already in model.")
        else:
            #TODO: Skip glucose (for one), has active transport from periplasm
            transport = Reaction(f"{met_id}-transport",
                                 lower_bound=0.0,
                                 upper_bound=abs(rate))
            transport._annotation["Experimental V_max (mM/hr/g)"] = abs(rate)
            transport.add_metabolites({met_ex : -1.0, met_c: 1.0})
            transport.gene_reaction_rule = gene_rule
            model.add_reactions([transport])

        # Constrain with fitted rate
        # TODO: Need to translate...?
        exchange._annotation["Experimental rate"] = rate

        # Assuming growth on glucose by default!
        if met_id != "Glucose[e]":
            exchange.bounds = (0, 0)
