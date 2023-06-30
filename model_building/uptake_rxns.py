import json
import re

from cobra.core import Metabolite

from utils.cobra_utils import change_compartment, get_or_create_exchange

UPTAKE_RATES_DATA = "parameters/uptake_rates/fitted_uptake_rates.json"
CARBON_SOURCE_IDS = "parameters/uptake_rates/carbon_source_ids.json"


def get_uptake_rates(uptake_rates_file=UPTAKE_RATES_DATA, carbon_sources_file=CARBON_SOURCE_IDS):
    with open(uptake_rates_file, "r") as f:
        uptake_rates = json.load(f)

    with open(carbon_sources_file, "r") as f:
        carbon_sources = json.load(f)

    result = {
        carbon_sources[metabolite]: rate
        for metabolite, rate in uptake_rates.items()
    }

    return result


def add_uptake_reactions(model, rates):
    # TODO: EX_glc reversibility is wrong given it should be an exchange reaction,
    # need to match other EX reactions?

    # Add uptake reactions one-by-one
    for met_id, rate in rates.items():
        # Check if exchange reaction already exists,
        # create exchange reaction if it does not
        try:
            exchange = get_or_create_exchange(model, met_id, verbose=True)
        except KeyError:
            print(f"{met_id} not found in model, failed to create exchange reaction.")
            continue

        # Constrain with fitted rate
        # TODO: Need to translate...?
        exchange._annotation["Experimental rate"] = rate

        # Assuming growth on glucose by default!
        if met_id != "Glucose[e]":
            exchange.bounds = (0, 0)
