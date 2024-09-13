import pickle
import re
import requests
import xmltodict
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from getpass import getpass
from cobra.io import read_sbml_model


RPOM_ORGID = "GCF_000011965"

def get_session(username=None, password=None):
    # Prompt for username and password if not provided
    if username is None:
        username = input("Username: ")
    if password is None:
        password = getpass("Password: ")
    
    # Get and return session
    s = requests.Session()
    r = s.post("https://websvc.biocyc.org/credentials/login/",
            data={"email": username, "password": password})
    r.raise_for_status()

    return s


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


def get_pathways(reaction, orgid=RPOM_ORGID, session=None):
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
    
    # Get pathways of reactions
    pathways = {}
    with tqdm(total=len(model.reactions)) as pbar:
        with ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(get_pathways, stems[reaction.id], RPOM_ORGID, (username, password)) : reaction
                for reaction in model.reactions}
            
            for future in as_completed(futures):
                pathways[futures[future].id] = future.result()
                pbar.update(1)

    # Save data
    with open("model_building/reaction_annotations.pkl", "wb") as f:
        pickle.dump({
            "stems" : stems,
            "pathways" : pathways
        }, f)



if __name__ == "__main__":
    main()
