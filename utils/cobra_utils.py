import re
from cobra.core import Metabolite, Reaction
from cobra.io import read_sbml_model
import networkx as nx
import heapq
import matplotlib.pyplot as plt
import numpy as np
from molmass import Formula


def widest_path(fluxes,
                source_met,
                target_met,
                ignore_ids={"WATER[c]",
                            "ATP[c]",
                            "ADP[c]",
                            "Pi[c]",
                            "PROTON[c]",
                            "NAD[c]",
                            "NADP[c]",
                            "NADH[c]",
                            "NADPH[c]", }):
    """Find the widest (least-constraining) path from source_met to target_met
    with the given fluxes. Paths through currency metabolites or certain reactions are blocked by specifying
    currency metabolites in the ignore_ids parameter.
    """
    def neighbors(node):
        if isinstance(node, Metabolite):
            # Return only reactions that consuming the metabolite
            result = set()
            for reaction in node.reactions:
                coeff = reaction.metabolites[node]
                if coeff * fluxes[reaction.id] < 0 and reaction.id not in ignore_ids:
                    result.add(reaction)
            return result
        elif isinstance(node, Reaction):
            # Return only metabolites producing the metabolite
            return set(met
                       for met, coeff in node.metabolites.items()
                       if coeff * fluxes[node.id] > 0
                       and met.id not in ignore_ids)
        else:
            raise ValueError(
                f"Expected a Metabolite or Reaction, got {type(node)}")

    def clean_path(path):
        # Utility function to remove extraneous nodes from the path found.
        # Traverse backwards from the target to the source in order to identify a single chain path
        node = target_met
        chain = [node]
        preds = list(path.predecessors(node))
        while preds:
            assert len(preds) == 1
            chain.extend(preds)
            node = preds[0]
            preds = list(path.predecessors(node))
        chain = list(reversed(chain))

        to_remove = set(path.nodes) - set(chain)
        path.remove_nodes_from(to_remove)
        return path

    # Resulting path will be given as a DiGraph
    path = nx.DiGraph()

    # priority queue: (-capacity, tiebreaker, (prev_node, node))
    _tiebreaker = 0
    pq = [(-float('inf'), _tiebreaker, (None, source_met))]
    visited = set()
    while pq:
        # Visit next node
        neg_capacity, _, (prev_node, node) = heapq.heappop(pq)

        # Ignore visited nodes
        if node in visited:
            continue
        visited.add(node)

        # Add to path
        path.add_node(node, node_type=type(node))
        if prev_node is not None:
            path.add_edge(prev_node, node, weight=-neg_capacity)

        # Stop if we've reached the target
        if node == target_met:
            return clean_path(path)

        # Move on to neighbors of node
        for neighbor in neighbors(node):
            # Take the capacity on this edge as the quantity per hour going through the reaction
            # (multiply flux * coefficient)
            # Ideally, we would convert this to mass, but some metabolites have generic formulas (e.g. with an R or X)
            met = neighbor if isinstance(neighbor, Metabolite) else node
            rxn = neighbor if isinstance(neighbor, Reaction) else node
            edge_capacity = abs(fluxes[rxn.id]
                                * rxn.metabolites[met])
            _tiebreaker += 1
            heapq.heappush(pq, (max(neg_capacity, -edge_capacity),
                           _tiebreaker, (node, neighbor)))

    # If we haven't returned by now, there is no path
    return None


def draw_path(path, mets_per_row = 4, jitter=2):
    pos = {}
    i, j = 0, 0
    direction = 1
    for node in nx.topological_sort(path):
        pos[node] = np.array([10*i,
                              -10*j + jitter*((i+j)%2)])
        
        if i + direction >= mets_per_row or i + direction < 0:
            
            j += 1
            direction *= -1
        else:
            i += direction
        

    nx.draw(path,
            pos=pos,
            labels={node: node.id if isinstance(node, Metabolite) else "" for node in path.nodes},
            with_labels=True)
    nx.draw_networkx_edge_labels(path, pos, edge_labels={(u, v): f"{path[u][v]['weight']:.2g}" for u, v in path.edges})


def reactions_of_metabolite(model, metID):
    return model.reactions.query(lambda r: metID in [met.id for met in r.metabolites],
                                 attribute=None)


def distance(model, o1, o2):
    queue = []
    if isinstance(o1, Metabolite):
        queue += [(r, 1) for r in reactions_of_metabolite(model, o1.id)]
    else:
        queue += [(m, 1) for m in o1.metabolites.keys()]
    
    seen = {o1: 0}
    while len(queue) > 0:
        next_item, steps = queue.pop()
        
        # Found the other item, return
        if next_item == o2:
            return steps
        
        # Do not get trapped in cycles
        if next_item in seen:
            continue
        seen[next_item] = steps

        # Add next items to queue
        if isinstance(next_item, Metabolite):
            queue += [(r, steps + 1) for r in reactions_of_metabolite(model, next_item.id)]
        else:
            queue += [(m, steps + 1) for m in next_item.metabolites.keys()]
    
    return float('inf')


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


def set_active_bound(reaction: Reaction, bound: float, abs_bounds=True) -> None:
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

    # TODO: This abs_bounds thing is really shoddy but I want to avoid changing set_active_bound everywhere it crops up for now
    if reaction.reactants:
        reaction.lower_bound = -abs(bound) if abs_bounds else bound
    elif reaction.products:
        reaction.upper_bound = abs(bound) if abs_bounds else bound


def set_fixed(reaction: Reaction, flux: float) -> Reaction:
    reaction.lower_bound = flux
    reaction.upper_bound = flux
    return reaction


def formula_str(reaction: Reaction) -> str:
    stoich = {
        metabolite : (metabolite.formula, coeff)
        for metabolite, coeff in reaction.metabolites.items()
    }

    arrow = "-->"
    if reaction.upper_bound == 0:
        if reaction.lower_bound == 0:
            arrow = "="
        else:
            arrow = "<--"
    elif reaction.reversibility:
        arrow = "<=>"

    left = " + ".join(f"{-coeff} {formula}" for _, (formula, coeff) in stoich.items() if coeff < 0)
    right = " + ".join(f"{coeff} {formula}" for _, (formula, coeff) in stoich.items() if coeff > 0)

    return f"{left} {arrow} {right}"
