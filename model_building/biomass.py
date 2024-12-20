import json

from cobra import Reaction, Metabolite
from cobra.io import load_model


MANUAL_METABOLITE_MATCHES = "model_building/manual_met_links.json"
HWA_BIOMASS = "model_building/reactions/hwa_biomass.json"


def get_ecoli_biomass():
    ecoli_model = load_model("iJO1366")

    # ID of biomass reaction - there may be a more future-proof way
    # of doing this (probably pulling from ecoli_model.objective), but fine for now...
    BIOMASS_RXN = "BIOMASS_Ec_iJO1366_core_53p95M"

    # Dictionary of metabolite : coefficient
    ecoli_biomass = ecoli_model.reactions.get_by_id(BIOMASS_RXN).metabolites

    return ecoli_biomass


def ecoli_biomass_to_rpom(r_pom_model, ecoli_biomass, manual_matches={}):
    # Prepare metabolite : coefficient dictionary
    rpom_biomass_rxn = {}

    # Add metabolites and coefficients from ecoli
    # into rpom_biomass_rxn one by one
    n_kegg_matches = 0
    n_manual_matches = 0
    for metabolite, coefficient in ecoli_biomass.items():
        # Get annotations to do linking
        kegg = metabolite.annotation.get("kegg.compound", "!!NO_KEGG_ID!!")
        formula = metabolite.formula
        compartment = metabolite.compartment

        # match metabolites by Kegg ID,
        # suggesting fall-backs based on name, formula
        if isinstance(kegg, list):
            # Sometimes, multiple kegg ids exist for metabolite in ecoli model -
            # need to query all, list unique matches found
            matching_metabolites = set()
            for id in kegg:
                query_result = r_pom_model.metabolites.query(
                    lambda m: m.annotation.get('Kegg ID') == id and
                    m.compartment == compartment
                )

                for m in query_result:
                    matching_metabolites.add(m)

            matching_metabolites = list(matching_metabolites)
        else:
            # kegg is a str, just need to search one id
            matching_metabolites = r_pom_model.metabolites.query(
                lambda m: m.annotation.get('Kegg ID') == kegg and
                m.compartment == compartment
            )

        # Add unique match to reaction if found
        if len(matching_metabolites) == 1:
            rpom_biomass_rxn[matching_metabolites[0]] = coefficient
            n_kegg_matches += 1
        # List matches if multiple found
        elif len(matching_metabolites) > 1:
            print()
            print("Multiple matches found for Kegg ID {kegg}.")
            for m in matching_metabolites:
                print(f"{metabolite.id} ({metabolite.name}) = {m.id} ({m.name})?")
        # No Kegg matches found, do fallback search by formula and list results
        else:
            print()
            print("No matches found in R. pom model for metabolite "
                  f"{metabolite.id} ({metabolite.name})"
                  f"\nwith Kegg ID {kegg}.")
            print("Falling back to search by formula...")

            matching_metabolites = r_pom_model.metabolites.query(
                lambda m: m.formula == formula and m.compartment == compartment
            )

            print(f"Found {len(matching_metabolites)} match(es).")

            for m in matching_metabolites:
                print(
                    f"- {metabolite.id} ({metabolite.name}) = {m.id} ({m.name})?")

        # Use manual match if supplied
        if metabolite.id in manual_matches:
            print(
                f"Using supplied manual match, {metabolite.id} = {manual_matches[metabolite.id]}.")

            matched_met = r_pom_model.metabolites.get_by_id(
                manual_matches[metabolite.id])
            rpom_biomass_rxn[matched_met] = coefficient
            n_manual_matches += 1

    print(
        f"\nSuccessfully linked {n_kegg_matches} metabolites by Kegg ID,"
        f"\nwith an additional {n_manual_matches} matches supplied manually"
        f"\nfor a total of {n_kegg_matches + n_manual_matches} / {len(ecoli_biomass)} successful links.\n")

    # Remove biotin - not produced by R pom
    biotin = r_pom_model.metabolites.get_by_id("BIOTIN[c]")
    del rpom_biomass_rxn[biotin]

    return rpom_biomass_rxn


def get_core_biomass_stoich():
    with open("data/objectives/core_biomass.json", "r") as f:
        core_biomass = json.load(f)
        biomass_objective = core_biomass["core_biomass"]

    return biomass_objective


def add_ecoli_core_biomass_to_model(model):
    biomass_stoich = get_core_biomass_stoich()
    biomass_stoich = {model.metabolites.get_by_id(
        met): coeff for met, coeff in biomass_stoich.items()}
    biomass_rxn = Reaction("RPOM_provisional_biomass",
                           "RPOM_provisional_biomass",
                           lower_bound=0,
                           upper_bound=1000)
    biomass_rxn.add_metabolites(biomass_stoich)

    # Add to model, set as objective
    model.add_reactions([biomass_rxn])
    model.objective = biomass_rxn


def add_ecoli_full_biomass_to_model(model):
    # Get ecoli biomass reaction stoichiometry
    ecoli_biomass = get_ecoli_biomass()

    # Load manual metabolite id links between E. coli and R. pom.
    with open(MANUAL_METABOLITE_MATCHES, "r") as f:
        manual_matches = json.load(f)

    # Get stoichiometry of biomass reaction
    biomass_stoich = ecoli_biomass_to_rpom(
        model, ecoli_biomass, manual_matches=manual_matches)

    # Create biomass reaction
    biomass_rxn = Reaction("RPOM_provisional_biomass",
                           "RPOM_provisional_biomass",
                           lower_bound=0,
                           upper_bound=1000)
    biomass_rxn.add_metabolites(biomass_stoich)

    # Add to model, set as objective
    model.add_reactions([biomass_rxn])
    model.objective = biomass_rxn


def add_hwa_biomass_to_model(model):
    with open(HWA_BIOMASS, "r") as f:
        biomass_reactions = json.load(f)

    # Create missing pseudo-metabolites
    biomass_stoich = biomass_reactions["Rpom_hwa_biomass"]
    for met in biomass_stoich:
        try:
            model.metabolites.get_by_id(met)
        except:
            comp = met[met.index("[")+1:met.index("]")]
            model.add_metabolites(Metabolite(met, compartment=comp))
            print(f"Added {met}.")

    # Create reactions
    for rxnid, stoich in biomass_reactions.items():
        reaction = Reaction(rxnid, rxnid)
        reaction.add_metabolites(
            {
                model.metabolites.get_by_id(metid): coeff
                for metid, coeff in stoich.items()
            }
        )
        model.add_reactions([reaction])
    
    model.objective = model.reactions.get_by_id("Rpom_hwa_biomass")
