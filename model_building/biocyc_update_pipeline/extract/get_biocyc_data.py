import pickle
import requests
import warnings
import xmltodict
from getpass import getpass
from utils.biocyc import get_session
from urllib.parse import quote_plus


def get_all_reactions(orgid, session=None, detail="full"):
    # Get all reactions for an organism from the BioCyc API
    if session is None:
        session = get_session()

    # Using ^^ instead of ^^^ includes reaction classes
    query = quote_plus(f"[x:x<-{orgid}^^reactions]")
    r = session.get(
        f"https://websvc.biocyc.org/xmlquery?{query}&detail={detail}")
    r.raise_for_status()
    return xmltodict.parse(r.text)["ptools-xml"]["Reaction"]


def wrap(slot_value):
    if not isinstance(slot_value, list):
        slot_value = [slot_value]
    return slot_value


def get_stoichiometry(reaction_data):
    # Utility function for retrieving stoichiometry from reaction data dictionary
    # (as obtained from biocyc API request)
    stoich = {}

    # Get reactants and products
    left = wrap(reaction_data.get("left", None))
    right = wrap(reaction_data.get("right", None))

    for side, sign in [(left, -1), (right, 1)]:
        for compound in side:
            if compound is None:
                continue

            # Get data for the compound, which may be one of several forms (based on my observations).
            # This is so unnecessarily complicated! I've seen the following forms:
            #
            # 1. A bare string, indicating a compound name, e.g. "a supercoiled duplex DNA"
            #    (probably not matching to any actual object in the database)
            #
            # 2. A dictionary with some of the following keys in some combination (together with example values):
            #        - "Compound"/"Protein"/"RNA" : {"@frameid": ..., ...}
            #               The most "standard" key, representing the schema class of the compound,
            #               and containing the frameid and other standard data.
            #
            #        - "#text" : "a deoxyribonucleic acid"
            #              Typically used to name compounds that are not actually objects in the database,
            #              e.g. electrons in redox half-reactions ("e<SUP>-</SUP>") or generic compounds ("a deoxyribonucleic acid", "an antibiotic", ...).
            #              Use this as a fallback for ID if the frameid is not present.
            #
            #        - "name-slot" : "N-Name" / "N-1-Name" / "N+1-Name" / ...,
            #              Used in polymerization reactions to indicate the length of the polymer.
            #              I noticed at least one instance of "N-2-NAME", with the "NAME" capitalized, so take note and
            #              consider case-insensitivity. I'll tack this on to the end of the compound id, whether from @frameid or #text.
            #
            #        - "coefficient" : {"@datatype": "integer", "#text": "3"} or {"@datatype": "string", "#text": "n-1"}
            #               The stoichiometric coefficient of the compound (always positive, so we need to add our own signs).
            #               May also be a string like "n" or "n-1", representing an undefined number of compounds
            #               (presumably in reaction classes for polymerization). We'll preserve such values, but recommend deletion down the line.
            #               If not present, defaults to 1.
            #
            #        - "compartment" : {"cco" : {..., "@frameid": ...}},
            #              The compartment of the compound, as a dictionary with the frameid of the compartment.
            #              Many of these will be the abstract compartments "CCO-IN" and "CCO-OUT", which really should be
            #              clarified through a "RXN-LOCATIONS" slot in the PGDB, though that slot does not appear to be accessible like this.
            #              We'll assume that CCO-IN is the cytoplasm and CCO-OUT is the periplasm (represented in our model as [c] and [p], respectively).
            #              The default is the cytoplasm ('CCO-CYTOSOL').

            compound_id = ""
            coefficient = sign  # Default coefficient is -1 for reactants, 1 for products
            if isinstance(compound, str):
                # Reactant is a compound name, so we'll use the name as the ID
                # and assume it's in the cytoplasm
                compound_id = f"{compound}[c]"
            else:
                # Reactant is a dictionary, so we'll try to extract the frameid and coefficient
                # First, look for a compound, protein, or RNA:
                found_id = False
                for key in ["Compound", "Protein", "RNA"]:
                    if key in compound:
                        compound_id = compound[key]["@frameid"]
                        found_id = True
                        break

                if not found_id:
                    # If we didn't find a frameid, use the #text field
                    compound_id = compound.get("#text", "")

                # For polymerization reactions, add the polymer length to the end of the compound id
                if "name-slot" in compound:
                    compound_id += f"__{compound['name-slot']}"

                # Get the compartment
                compartment = compound.get("compartment", {}).get(
                    "cco", {}).get("@frameid", "CCO-CYTOSOL")
                if compartment == "CCO-CYTOSOL":
                    compartment = "c"
                elif compartment == "CCO-IN":
                    compartment = "c"
                elif compartment == "CCO-OUT":
                    compartment = "p"
                else:
                    warnings.warn(
                        f"Unknown compartment frameid: {compartment}")

                compound_id += f"[{compartment}]"

                # Get the stoichiometric coefficient (defaults to -1, see above)
                # and make it negative to indicate a reactant
                if "coefficient" in compound:
                    if compound["coefficient"]["@datatype"] == "integer":
                        coefficient = sign * \
                            int(compound["coefficient"]["#text"])
                    else:
                        coefficient = sign * compound["coefficient"]["#text"]

            stoich[compound_id] = coefficient

    return stoich


def get_metabolites_from(orgid, met_ids, session=None, detail="full"):
    # Get all metabolites in met_ids from specified organism using BioCyc API
    if session is None:
        session = get_session()
    
    metabolites = []
    for met_id in met_ids:
        try:
            r = session.get(
                f"https://websvc.biocyc.org/getxml?id={orgid}:{quote_plus(met_id)}&detail={detail}")
            r.raise_for_status()
            met_data = xmltodict.parse(r.text)["ptools-xml"]
            for type_key in ["Compound", "Protein", "RNA"]:
                if type_key in met_data:
                    met_data = met_data[type_key]
                    break
            metabolites.append(met_data)
        except requests.HTTPError as e:
            print(f"Error retrieving metabolite {met_id}: {e}")
            continue
        except KeyError as e:
            print(f"Error parsing metabolite {met_id}: {list(xmltodict.parse(r.text)["ptools-xml"].keys())}")
            continue

    return metabolites


def get_all_proteins(orgid, session=None, detail="full"):
    # Get all proteins for an organism from the BioCyc API
    if session is None:
        session = get_session()

    # Using ^^^ instead of ^^ excludes protein classes
    query = quote_plus(f"[x:x<-{orgid}^^^proteins]")
    r = session.get(
        f"https://websvc.biocyc.org/xmlquery?{query}&detail={detail}")
    r.raise_for_status()
    return xmltodict.parse(r.text)["ptools-xml"]["Protein"]


def get_all_genes(orgid, session=None, detail="full"):
    # Get all genes for an organism from the BioCyc API
    if session is None:
        session = get_session()

    # Using ^^^ instead of ^^ excludes gene classes
    query = quote_plus(f"[x:x<-{orgid}^^^genes]")
    r = session.get(
        f"https://websvc.biocyc.org/xmlquery?{query}&detail={detail}")
    r.raise_for_status()
    return xmltodict.parse(r.text)["ptools-xml"]["Gene"]


def main():
    DB1 = "RUEGERIA_POMEROYI_DSS3"
    DB2 = "GCF_000011965"

    s = get_session()

    # Get reaction data from databases, save off
    for database in [DB1, DB2]:
        print(f"Getting reactions for {database}...")
        reactions = get_all_reactions(database, s)
        with open(f"model_building/biocyc_update_pipeline/data/raw_reaction_data__{database}.pkl", "wb") as f:
            pickle.dump(reactions, f)
    
        # Create a list of all metabolites in the reactions
        all_metabolites = set()
        for reaction in reactions:
            stoich = get_stoichiometry(reaction)
            all_metabolites.update([m.split('[')[0].split("__N")[0] for m in stoich.keys()])
    
        # Get metabolite data from database, save off
        print(f"Getting metabolites for {database}...")
        metabolites = get_metabolites_from(database, list(all_metabolites), s)

        with open(f"model_building/biocyc_update_pipeline/data/raw_metabolite_data__{database}.pkl", "wb") as f:
            pickle.dump(metabolites, f)

        # Get protein data from database, save off
        print(f"Getting proteins for {database}...")
        proteins = get_all_proteins(database, s)
        with open(f"model_building/biocyc_update_pipeline/data/raw_protein_data__{database}.pkl", "wb") as f:
            pickle.dump(proteins, f)

        # Get gene data from database, save off
        print(f"Getting genes for {database}...")
        genes = get_all_genes(database, s)
        with open(f"model_building/biocyc_update_pipeline/data/raw_gene_data__{database}.pkl", "wb") as f:
            pickle.dump(genes, f)


if __name__ == "__main__":
    main()
