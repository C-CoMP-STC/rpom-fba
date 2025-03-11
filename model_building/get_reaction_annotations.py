from collections import defaultdict
import pickle
import re
import requests
import xmltodict
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from getpass import getpass
from cobra.io import read_sbml_model
from urllib.parse import quote_plus
from utils.biocyc import get_session


RPOM_ORGID = "GCF_000011965"


def get_stem(reaction_full_id):
    # Only consider "specified" reactions
    if "//" not in reaction_full_id:
        return reaction_full_id
    
    # Okay so, these ids are really badly constructed and not
    # consistent at all. So we're going to have to do some
    # string manipulation to get the stem

    # All "bad" reactions have one of the following substrings:
    # 'RXN1G48': 321, 'RXN': 160, 'RXN0': 9, 'RXN66': 3, 'RXN1RFH': 2

    # Reactions with RXN1G48 are probably not in Biocyc (seem to be transport reactions)
    # Anyway, they're formed like TRANS-RXN1G48-##-<specification>
    if "RXN1G48" in reaction_full_id:
        return re.match(r"(^.*RXN1G48-\d+)-.*$", reaction_full_id).group(1)
    
    # Reactions with RXN0 are formed like ...RXN0-####-<specification>
    if "RXN0" in reaction_full_id:
        return re.match(r"(^.*RXN0-\d+)-.*$", reaction_full_id).group(1)
    
    # Reactions with RXN66 are formed like ...RXN66-####-<specification>
    if "RXN66" in reaction_full_id:
        return re.match(r"(^.*RXN66-\d+)-.*$", reaction_full_id).group(1)
    
    # Reactions with RXN1RFH are formed like ...RXN1RFH-####-<specification>
    if "RXN1RFH" in reaction_full_id:
        return re.match(r"(^.*RXN1RFH-\d+)-.*$", reaction_full_id).group(1)
    
    # Reactions with RXN are (apparently?) formed like ...RXN-<specification>
    if "RXN" in reaction_full_id:
        return re.match(r"(^.*RXN)-.*$", reaction_full_id).group(1)


def genes_of_reaction(reaction, orgid=RPOM_ORGID, session=None):
    # Use session if provided
    # (can be Session object or tuple of username, password),
    # or create one if not
    if isinstance(session, requests.Session):
        s = session
    elif isinstance(session, tuple):
        username, password = session
        s = get_session(username, password)
    else:
        s = get_session()
    
    # Request genes of reaction
    r = s.get(f"https://websvc.biocyc.org/apixml?fn=genes-of-reaction&id={orgid}:{reaction}&detail=none")
    if r.status_code == 404:
        return []
    
    # Clean up response
    genes = xmltodict.parse(r.text)["ptools-xml"].get("Gene", [])
    if isinstance(genes, dict):
        genes = [genes]
    genes = [gene["@frameid"] for gene in genes]

    return genes


def get_pathway_data(orgid=RPOM_ORGID, session=None):
    # Use session if provided
    # (can be Session object or tuple of username, password),
    # or create one if not
    if isinstance(session, requests.Session):
        s = session
    elif isinstance(session, tuple):
        username, password = session
        s = get_session(username, password)
    else:
        s = get_session()

    # Using ^^^ instead of ^^ to get only instances (not classes)
    query = quote_plus(f"[x:x<-{orgid}^^^pathways]")
    r = s.get(f"https://websvc.biocyc.org/xmlquery?{query}")
    r.raise_for_status()

    # Clean up response
    pathways = xmltodict.parse(r.text)["ptools-xml"].get("Pathway", [])

    # Key by frameid
    pathway_data = {
        pwy["@frameid"] : pwy
        for pwy in pathways
    }

    # Function to recursively get all reactions in a pathway
    def get_reactions(pwy_dict):
        reaction_list = pwy_dict.get("reaction-list", {})
        result = []

        # If there are subpathways, get reactions from them recursively
        if "Pathway" in reaction_list:
            subpathways = reaction_list["Pathway"]
            if not isinstance(subpathways, list):
                subpathways = [subpathways]
            
            for subpwy in subpathways:
                subpwy_id = subpwy["@frameid"]
                result.extend(get_reactions(pathway_data[subpwy_id]))
        
        # Get reactions from this pathway
        reactions = reaction_list.get("Reaction", [])
        if not isinstance(reactions, list):
            reactions = [reactions]
        for reaction in reactions:
            result.append(reaction["@frameid"])

        return result
    
    # Get all reactions in each pathway
    pathway_data = {
        pwy_id : {
            "common-name": pwy_dict.get("common-name", {}).get("#text", ""),
            "reactions": get_reactions(pwy_dict)
        }
        for pwy_id, pwy_dict in pathway_data.items()
    }

    return pathway_data


def get_pathways_of_genes_of_reaction(reaction, orgid=RPOM_ORGID, session=None):
    """DEPRECATED: 
    Unfortunately, there is no direct way to get the pathways of a reaction
    with BioCyc API functions. This method is a workaround that gets the genes
    of a reaction, then gets the pathways of each gene. However, this approach
    is too inclusive, since a gene can catalyze multiple reactions belonging to
    different pathways. This method is kept here for reference, but is not used.

    Instead, use the get_pathway_data function to get all pathways in the organism,
    then reactions in each pathway directly from that data.
    """
    # Use session if provided
    # (can be Session object or tuple of username, password),
    # or create one if not
    if isinstance(session, requests.Session):
        s = session
    elif isinstance(session, tuple):
        username, password = session
        s = get_session(username, password)
    else:
        s = get_session()

    s = get_session(username, password)
    genes = genes_of_reaction(reaction, orgid, s)

    
    pathway_ids = set()
    for gene in genes:
        r = s.get(f"https://websvc.biocyc.org/apixml?fn=pathways-of-gene&id={orgid}:{gene}&detail=low")
        gene_pathways = xmltodict.parse(r.text)["ptools-xml"].get("Pathway", [])
        if isinstance(gene_pathways, dict):
            gene_pathways = [gene_pathways]
        pathway_ids.update([path["@frameid"] for path in gene_pathways])
    
    return set(pathway_ids)


def main():
    MODEL = "model/Rpom_05.xml"
    model = read_sbml_model(MODEL)

    # Get bicyc credentials
    username = input("Username: ")
    password = getpass("\nPassword: ")

    # Check login works
    s = requests.Session()
    r = s.post("https://websvc.biocyc.org/credentials/login/",
            data={"email": username, "password": password})
    r.raise_for_status()

    print("Successfully logged in.\n")

    # Get stems of reactions
    stems = {
        reaction.id: get_stem(reaction.id) for reaction in model.reactions
    }

    # Get all pathways in R. pom
    pathway_data = get_pathway_data(RPOM_ORGID, (username, password))

    # Get pathways of reactions, as a {reaction: [pathways]} dict.
    # Note that in order to account for instance reactions, we use the *stem*
    # to check for membership in pathway.
    pathways = defaultdict(list)
    # remaining = set(model.reactions)
    # for pwy_id, pwy_data in tqdm(pathway_data.items()):
    #     for reaction in pwy_data["reactions"]:
    #         pathways[reaction].append(pwy_id)
    #         remaining.discard(reaction)
    
    # # Remaining reactions are not in any pathway
    # for reaction in remaining:
    #     pathways[reaction] = []
    for reaction, stem in tqdm(stems.items()):
        for pwy_id, pwy_data in pathway_data.items():
            if stem in pwy_data["reactions"]:
                pathways[reaction].append(pwy_id)

    # Save pathway data
    with open("model_building/pathway_data.pkl", "wb") as f:
        pickle.dump(pathway_data, f)

    # Save reaction annotations
    with open("model_building/reaction_annotations.pkl", "wb") as f:
        pickle.dump({
            "stems" : stems,
            "pathways" : pathways
        }, f)



if __name__ == "__main__":
    main()
