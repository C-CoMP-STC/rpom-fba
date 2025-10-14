from getpass import getpass
import os
import pandas as pd
import numpy as np

from cobra.core import Reaction, Model
from cobra.io import read_sbml_model
import requests
from tqdm import tqdm

from utils.cobra_utils import get_or_create_exchange


def to_paintable(model : Model, collection : str = "reactions", pbar : bool = False, **data):
    """Converts the specified collection within a COBRA model to a dataframe,
    which can be painted onto the model in BioCyc Cellular Overview.

    Args:
        model (Model): The COBRA model.
        collection (str, optional): The collection of objects within the model to include.\
            Can be "reactions", "metabolites", or "genes". Defaults to "reactions".
        pbar (bool, optional): Whether to display a progress bar. Defaults to False.
        **data: Data items to include in the dataframe.\
            Each key becomes the name of a column in the dataframe, while the items become the column content.\
            Each item can be a list of values, a single value, or a callable that takes an object from the collection\
            and returns a value.

    Returns:
        DataFrame : A dataframe with the specified collection and data items.
    
    Raises:
        ValueError: If an invalid collection is specified.
    """

    # Get appropriate collection:
    match collection:
        case "reactions":
            objects = model.reactions
        case "metabolites":
            objects = model.metabolites
        case "genes":
            objects = model.genes
        case _:
            raise ValueError(f"Invalid collection: {collection}")

    # For data items that are callable, call them on each element of the collection:
    for colname, value in data.items():
        if callable(value):
            data[colname] = [value(item) for item in objects]

    # For each item in the collection, add a row with the item id
    # and corresponding data entries:
    result = []

    for i, item in enumerate(objects if not pbar else tqdm(objects)):
        result.append({
            "id": item.id if not isinstance(item, Reaction) else item.notes["stem"],
            **{
                colname: value[i] if isinstance(value, list) else value
                for colname, value in data.items()
            }
        })
    
    return pd.DataFrame(result)


def get_session(username=None, password=None, n_tries=3):
    # Prompt for username and password if not provided
    if username is None:
        username = input("Username: ")
    if password is None:
        password = getpass("Password: ")
    
    # Get and return session
    for _ in range(n_tries):
        try:
            s = requests.Session()
            r = s.post("https://websvc.biocyc.org/credentials/login/",
                    data={"email": username, "password": password})
            r.raise_for_status()
            return s
        except requests.HTTPError as e:
            print(f"Error: {e}")
            username = input("Username: ")
            password = getpass("Password: ")

    return s


def main():
    # Load model and retrieve exchange reactions
    model = read_sbml_model("model/Rpom_05.xml")
    ex_glc = model.reactions.get_by_id("EX_glc")
    ex_ace = get_or_create_exchange(model, "ACET[e]")

    # Turn on maintenance
    atpm = model.reactions.get_by_id("ATPM")
    atpm.bounds = (25, 25)

    # Get fluxes for growth on glucose
    with model:
        ex_glc.lower_bound = -3
        sol = model.optimize()
        flux_glc = sol.fluxes / sol.objective_value

    # Get fluxes for growth on acetate
    with model:
        ex_ace.lower_bound = -9
        sol = model.optimize()
        flux_ace = sol.fluxes / sol.objective_value

    rxn_fc = flux_ace / flux_glc

    # Create output directory
    os.makedirs("out/biocyc", exist_ok=True)

    # Generate paintable for reaction presence
    presence_df = to_paintable(model,
                               collection="reactions",
                               pbar=True,
                               presence=1)
    presence_df.to_csv("out/biocyc/reaction_presence.tsv", index=False, sep="\t")

    # Generate paintable for reaction fluxes
    flux_df = to_paintable(
        model,
        collection="reactions",
        pbar=True,
        flux_glc = lambda r: flux_glc.get(r.id, 0),
        flux_ace = lambda r: flux_ace.get(r.id, 0)
    )
    flux_df = flux_df[(flux_df["flux_glc"] != 0) | (flux_df["flux_ace"] != 0)]
    flux_df.to_csv("out/biocyc/reaction_flux.tsv", index=False, sep="\t")

    # Generate paintable for reaction fold changes
    fc_df = to_paintable(
        model,
        collection="reactions",
        pbar=True,
        rxn_fc = lambda r: rxn_fc.get(r.id, 0)
    )
    fc_df = fc_df[fc_df["rxn_fc"] != 0]
    fc_df.to_csv("out/biocyc/reaction_fc.tsv", index=False, sep="\t")


    # flux_glc_df = to_paintable(model,
    #                            collection="reactions",
    #                            pbar=True,
    #                            flux_glc=lambda r: flux_glc.get(r.id, 0))
    # flux_glc_df = flux_glc_df[flux_glc_df["flux_glc"] != 0]
    # flux_glc_df.to_csv("out/biocyc/reaction_flux_glc.tsv", index=False, sep="\t")

    # # Generate paintable for reaction fluxes on acetate
    # flux_ace_df = to_paintable(model,
    #                            collection="reactions",
    #                            pbar=True,
    #                            flux_ace=lambda r: flux_ace.get(r.id, 0))
    # flux_ace_df = flux_ace_df[flux_ace_df["flux_ace"] != 0]
    # flux_ace_df.to_csv("out/biocyc/reaction_flux_ace.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()
