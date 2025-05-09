import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from itertools import count, chain


def formula_color(elements):
    element_colors = {
        "C" : np.array([0, 0, 1]),
        "N" : np.array([0, 1, 0]),
        "P" : np.array([1, 0, 1]),
        "H" : np.array([1, 0, 0]),
        "O" : np.array([0, 1, 0.5])
    }

    denom = sum(elements.values())

    result = sum(
        (elements.get(elem, 0) / denom) * c
        for elem, c in element_colors.items()
    )
    if any(result > 1):
        result /= result.sum()
    return result


def uptakes_and_secretions(model, boundary_fluxes, sol, initial_biomass, CUE_VOLUME, initial_glucose, initial_acetate):
    fig, ax = plt.subplots()

    sum_secreting = 0
    sum_uptake = 0
    for _, (rxn, flux, secreting, molmass) in boundary_fluxes.iterrows():
        flux = abs(flux) * molmass  # flux is in mmol/hr, so multiplying by molmass results in mg/hr
        if flux == 0:
            continue

        ax.bar(int(secreting),
            flux,
            width=0.8,
            bottom=sum_secreting if secreting else sum_uptake,
            align="center",
            color=formula_color(list(model.reactions.get_by_id(rxn).metabolites.keys())[0].elements),
            label=str(rxn))

        if secreting:
            sum_secreting += flux
        else:
            sum_uptake += flux

    ax.bar(1,
            sol.objective_value * initial_biomass * CUE_VOLUME.to("L").magnitude * 1000,
            width=0.8,
            bottom=sum_secreting,
            color="gray",
            align="center",
            label="biomass")
    # ax.hlines(0, 0, 1, ["k"])
    ax.set_xticks([0, 1], ["uptake", "secretions"])
    ax.set_ylabel("flux (mg / hr)")
    ax.set_title(f"{initial_glucose} mM Glucose, {initial_acetate} mM Acetate")
    fig.legend()

    return fig, ax


def flux_color(flux, c0 = (0.75, 0.75, 0.75), c_minus = (0, 0, 1), c_plus = (1, 0, 0)):
    c_end = c_plus
    if flux < 0:
        flux = -flux
        c_end = c_minus
    theta = (10 * flux) / (10 * flux + 1)
    return theta * np.array(c_end) + (1 - theta) * np.array(c0)


def plot_pathway(model, metabolite_graph, reaction_list, ax=None, sol=None):
    if ax is None:
        _, ax = plt.subplots()

    def rxn_flux(rxn):
        return rxn.flux if sol is None else sol.fluxes[rxn.id]

    for rxn_id in reaction_list:
        # Add node for reaction
        rxn = model.reactions.get_by_id(rxn_id)
        metabolite_graph.add_node(rxn_id,
                            reaction=True,
                            flux=rxn_flux(rxn))

        # Add edges to metabolites and position node
        met_positions = []
        for met, coeff in rxn.metabolites.items():
            if met.id in metabolite_graph.nodes:
                if coeff > 0:
                    metabolite_graph.add_edge(rxn_id, met.id, rxn=rxn_id, flux=rxn_flux(rxn))
                else:
                    metabolite_graph.add_edge(met.id, rxn_id, rxn=rxn_id, flux=rxn_flux(rxn))
                met_positions.append(metabolite_graph.nodes[met.id]["pos"])
        metabolite_graph.nodes[rxn_id]["pos"] = np.array(met_positions).mean(axis=0)

    nx.draw(metabolite_graph,
            with_labels=True,
            pos={n: metabolite_graph.nodes[n]["pos"] for n in metabolite_graph.nodes},
            node_size=[0 if metabolite_graph.nodes[n].get("reaction", False) else 100
                    for n in metabolite_graph.nodes],
            edge_color=[flux_color(metabolite_graph.edges[e]["flux"])
                        for e in metabolite_graph.edges],
            labels={n: f"{metabolite_graph.nodes[n]['flux']:.2e}" if metabolite_graph.nodes[n].get("reaction", False) else n
                for n in metabolite_graph.nodes},
            width=5,
            ax=ax)

    return ax


def plot_metabolite_fluxes(model, metabolite_id, ax=None, sol=None, label_reactions="fluxes", include_zeros=True, top_N=None):
    if ax is None:
        _, ax = plt.subplots()

    if top_N is None:
        top_N = len(model.reactions)

    # Initialize graph
    g = nx.DiGraph()
    metabolite = model.metabolites.get_by_id(metabolite_id)
    node_id = count()
    id_0 = next(node_id)
    g.add_node(id_0, label=metabolite_id, pos=(0, 0))

    # Put metabolite reactions in sorted order, only take up to top_N
    reactions = sorted(metabolite.reactions, key=lambda rxn: abs(rxn.flux), reverse=True)[:top_N]

    # Add reactions, metabolites to graph
    input_metabolites = []
    output_metabolites = []
    for reaction in reactions:
        flux = reaction.flux if sol is None else sol.fluxes[reaction.id]

        if not include_zeros and flux == 0:
            continue

        coeff = reaction.metabolites[metabolite]
        producing = flux * coeff > 0

        rxn_id = next(node_id)
        g.add_node(rxn_id,
                   reaction=True,
                   label=f"{flux:.2g}" if label_reactions == "fluxes" else reaction.id,
                   name=reaction.id)

        if producing:
            g.add_edge(rxn_id, id_0, flux = flux)
            for met, coeff in reaction.metabolites.items():
                # Skip central metabolite and anything else produced
                if met == metabolite or flux * coeff > 0:
                    continue
                met_id = next(node_id)
                input_metabolites.append(met_id)
                g.add_node(met_id, label=met.id, align="right")
                g.add_edge(met_id, rxn_id, flux = flux)

        else:
            g.add_edge(id_0, rxn_id, flux = flux)
            for met, coeff in reaction.metabolites.items():
                # Skip central metabolite and anything else consumed
                if met == metabolite or flux * coeff < 0:
                    continue
                met_id = next(node_id)
                output_metabolites.append(met_id)
                g.add_node(met_id, label=met.id, align="left")
                g.add_edge(rxn_id, met_id, flux = flux)

    # Position metabolite nodes
    n_inputs = len(input_metabolites)
    n_outputs = len(output_metabolites)
    for i, met in enumerate(input_metabolites):
        g.nodes[met]["pos"] = (-20, 10 * i - 5 * n_inputs)

    for i, met in enumerate(output_metabolites):
        g.nodes[met]["pos"] = (20, 10 * i - 5 * n_outputs)

    # Position reaction nodes
    for n in g.nodes:
        if g.nodes[n].get("reaction", False):
            # weighted average of metabolite positions, weighting (0,0) more heavily
            # (all mets get weight 1 except (0,0), with weight n - 1)
            neighbor_pos = np.array([g.nodes[m]["pos"]
                                          for m in chain(g.predecessors(n), g.successors(n))]
                                         )
            g.nodes[n]["pos"] = neighbor_pos.sum(axis=0) / (2 * neighbor_pos.shape[0] - 2)

    nx.draw_networkx(g,
                     pos={n: g.nodes[n]["pos"] for n in g.nodes},
                     labels={
                            n: g.nodes[n].get("label", "")
                            for n in g.nodes
                         },
                     #  horizontalalignment={n: g.nodes[n].get("align", "center") for n in g.nodes},
                     node_size = [0 if g.nodes[n].get("reaction", False) else 100 for n in g.nodes],
                     edge_color=[flux_color(g.edges[e]["flux"]) for e in g.edges],
                     width=2,
                     with_labels=True,
                     font_size=8,
                     ax=ax)

    return g, ax


def main():
    pass

if __name__ == "__main__":
    main()
