import re
from cobra.core import Metabolite, Reaction
from cobra.io import read_sbml_model


def reactions_of_metabolite(model, metID):
    return model.reactions.query(lambda r: metID in [met.id for met in r.metabolites],
                                 attribute=None)


def change_compartment(metabolite_id, new_compartment):
    return re.sub(r'\[.*?\]', f'[{new_compartment}]', metabolite_id)


def get_or_create_external_metabolite(model, cytoplasmic_metabolite_id):
    external_id = change_compartment(cytoplasmic_metabolite_id, "e")

    # Get cytoplasmic version of metabolite
    # (presumed to exist, will raise a KeyError otherwise)
    cyto_met = model.metabolites.get_by_id(cytoplasmic_metabolite_id)

    try:
        return model.metabolites.get_by_id(external_id)
    except KeyError:
        external_met = Metabolite(
            external_id,
            cyto_met.formula,
            cyto_met.name,
            cyto_met.charge,
            "e"
        )
        model.add_metabolites([external_met])

        return external_met


def get_or_create_exchange(model, external_metabolite_id, verbose=False):
    metabolite = get_or_create_external_metabolite(
        model, change_compartment(external_metabolite_id, "c"))

    # Check if exchange reaction already exists,
    # create exchange reaction if it does not
    exchange = [rxn for rxn in model.exchanges if metabolite in rxn.reactants]

    if len(exchange) > 1:
        raise ValueError(
            f"More than one exchange found for {external_metabolite_id}!")
    elif len(exchange) == 1:
        exchange = exchange[0]
    else:
        if verbose:
            print(
                f"Exchange reaction for {external_metabolite_id} not found, creating exchange reaction")
        exchange = model.add_boundary(metabolite, type="exchange")

    return exchange


def is_producing(reaction: Reaction) -> bool:
    """Determine if boundary reaction permits flux towards creating metabolites.

    Parameters
    ----------
    reaction: cobra.Reaction

    Returns
    -------
    bool
        True if reaction produces metaoblites and has upper_bound above 0
        or if reaction consumes metabolites and has lower_bound below 0 (so
        could be reversed).
    """
    return (bool(reaction.products) and (reaction.upper_bound > 0)) or (
        bool(reaction.reactants) and (reaction.lower_bound < 0)
    )


def get_active_bound(reaction: Reaction) -> float:
    """For an active boundary reaction, return the relevant bound.

    Parameters
    ----------
    reaction: cobra.Reaction

    Returns
    -------
    float:
        upper or minus lower bound, depenending if the reaction produces or
        consumes metaoblties.
    """
    if reaction.reactants:
        return -reaction.lower_bound
    elif reaction.products:
        return reaction.upper_bound


def set_active_bound(reaction: Reaction, bound: float) -> None:
    """Set active bound for boundary reaction (i.e., the bound constraining
    production of metabolites).

    Parameters
    ----------
    reaction: cobra.Reaction
        Reaction to set
    bound: float
        Value to set bound to. The bound is reversed and set as lower bound
        if reaction has reactants (metabolites that are consumed). If reaction
        has reactants, it seems the upper bound won't be set.
    """
    if reaction.reactants:
        reaction.lower_bound = -bound
    elif reaction.products:
        reaction.upper_bound = bound


def set_fixed(reaction: Reaction, flux: float) -> Reaction:
    reaction.lower_bound = flux
    reaction.upper_bound = flux
    return reaction


def path_to(model, met_from, met_to):
    pass


def diff_models(model1, model2):
    pass


def main():
    model1 = read_sbml_model("clean_models/Rpom_0.xml")
    model2 = read_sbml_model("clean_models/Rpom_02.xml")

    diff_models(model1, model2)


if __name__ == "__main__":
    main()
