from cobra.io import load_model

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

    return rpom_biomass_rxn
