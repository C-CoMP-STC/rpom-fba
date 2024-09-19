import pandas as pd

from cobra.core import Reaction, Model
from cobra.io import read_sbml_model
from tqdm import tqdm

from cobra_utils import get_or_create_exchange


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
            "id": item.id if not isinstance(item, Reaction) else item.annotation["stem"],
            **{
                colname: value[i] if isinstance(value, list) else value
                for colname, value in data.items()
            }
        })
    
    return pd.DataFrame(result)


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
        ex_glc.lower_bound = -10
        sol = model.optimize()
        flux_glc = sol.fluxes


    # Get fluxes for growth on acetate
    with model:
        ex_ace.lower_bound = -10
        sol = model.optimize()
        flux_ace = sol.fluxes


    # Generate paintable
    df = to_paintable(model,
                      collection = "reactions",
                      pbar=True,
                      presence = 1,
                      flux_glc = lambda r: flux_glc.get(r.id, 0),
                      flux_ace = lambda r: flux_ace.get(r.id, 0)
    )

    df.to_csv("out/biocyc/reaction_paintable.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()
