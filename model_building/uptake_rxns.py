import json

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
    # need to match other EX reactions

    # Clear existing exchange reactions?

    # Add uptake reactions one-by-one
    for met_id, rate in rates.items():

        # Check that metabolite is already in model
        try:
            metabolite = model.metabolites.get_by_id(met_id)
        except KeyError:
            # TODO: Fix missing metabolites
            print(f"{met_id} not found in model")
            continue

        # Check if exchange reaction already exists,
        # create exchange reaction if it does not
        exchange = model.exchanges.query(lambda r: metabolite in r.reactants)
        if len(exchange) == 0:
            print(f"Exchange reaction for {met_id} not found, creating exchange reaction")
            exchange = model.add_boundary(metabolite, type="exchange")

        # Constrain with fitted rate
        # Need to translate

        exchange.lower_bound = rate
