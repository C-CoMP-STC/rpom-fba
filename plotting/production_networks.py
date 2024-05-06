import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import networkx as nx
from collections import defaultdict


def active_graph(model):
    """Returns a reaction-metabolite bigraph for the given model, including only reactions with nonzero flux
    and metabolites connected to such reactions.

    Args:
        model (cobra.core.model.Model): the model for which to construct the bigraph.
    Returns:
        g (networkx.DiGraph): the bigraph.
    """
    g = nx.DiGraph()
    for reaction in model.reactions:
        if reaction.flux == 0:
            continue

        g.add_node(reaction.id, nodetype="reaction",
                   reaction_str=reaction.reaction, flux=abs(reaction.flux))
        for met, coeff in reaction.metabolites.items():
            g.add_node(met.id, nodetype="metabolite")
            flux_towards = coeff * reaction.flux
            n_from, n_to = (met.id, reaction.id) if flux_towards < 0 else (
                reaction.id, met.id)
            g.add_edge(n_from, n_to, flux=abs(reaction.flux))

    return g


def highest_flux_path(active_graph, source, target, INF_FLUX=1e3):
    return nx.shortest_path(active_graph, source, target, weight=lambda x, y, d: INF_FLUX-d["flux"])


def production_paths(active_graph, sources, targets, INF_FLUX=1e3):
    all_nodes = set()
    for s in sources:
        for t in targets:
            try:
                all_nodes.update(
                    highest_flux_path(active_graph, s, t, INF_FLUX))
            except nx.NetworkXException:
                continue
    return active_graph.subgraph(all_nodes)


def layout_production_subgraph(production_subgraph, inputs):
    inputs = [n for n in inputs if n in production_subgraph.nodes]

    # Assign a "depth" to each node based on how far it is from an input
    for depth, nodes in enumerate(nx.bfs_layers(production_subgraph, inputs)):
        nx.set_node_attributes(production_subgraph, {
                               n: depth for n in nodes}, "depth")
        for i, node in enumerate(nodes):
            production_subgraph.nodes[node]["pos"] = (15 * i, depth * 10)
    
    return production_subgraph


def plot_production_paths(production_subgraph, inputs, label=None, ax=None):
    if ax is None:
        _, ax = plt.subplots()

    layout_production_subgraph(production_subgraph, inputs)

    # If label is None, use node name, else if Node is a function(node_id, data), call the function
    if label is None:
        labels = {
            n: n if production_subgraph.nodes[n]["nodetype"] == "metabolite" else n[:10]
            for n in production_subgraph.nodes
        }
    elif callable(label):
        labels = {
            n: label(n, production_subgraph.nodes[n])
            for n in production_subgraph.nodes
        }
    else:
        raise ValueError(f"Unrecognized value for label, {label}")

    def color(flux): return 1 - (1 / (1 + np.exp(-flux)))
    nx.draw_networkx(production_subgraph,
                     pos=nx.get_node_attributes(production_subgraph, "pos"),
                     labels=labels,
                     node_size=[0 if production_subgraph.nodes[n]["nodetype"] == "reaction"
                                else 100
                                for n in production_subgraph.nodes],
                     edge_color=[[color(production_subgraph.edges[e]["flux"])]
                                 * 3 for e in production_subgraph.edges],
                     ax=ax
                     )

    return ax


def plot_network_as_plotly(G):
    pos = nx.get_node_attributes(G, "pos")

    node_x = [x for _, (x, y) in pos.items()]
    node_y = [y for _, (x, y) in pos.items()]
    nodetype = list(nx.get_node_attributes(G, "nodetype").values())

    edge_trace = []

    def color(flux):
        shade = 255 * (1 - (1 / (1 + np.exp(-flux))))
        return f"rgb({shade}, {shade}, {shade})"

    # Create edges as lines
    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]

        edge_trace.append(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None],
            line=dict(width=2, color=color(data["flux"])),
            mode='lines',
            hoverinfo='none'
        ))

    # Create nodes as a scatter trace
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=[node if data["nodetype"] == "metabolite" else ""
              for node, data in G.nodes(data=True)],
        textposition="top center",
        customdata=[
            node
            if data["nodetype"] == "metabolite"
            else f"{node}<br>{pathways_of_reaction(node)}"
            for node, data in G.nodes(data=True)],
        hovertemplate="%{customdata}",
        marker=dict(
            showscale=False,
            color='#1f78b4',
            size=[10 if nt == "metabolite" else 5 for nt in nodetype],
            line_width=2))

    # Create the figure to be displayed
    fig = go.Figure(data=[*edge_trace, node_trace],
                    layout=go.Layout(
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False,
                                   showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    # Display the figure
    return fig

def pathways_of_reaction(reaction_id, pathways_data="data/All-reactions-of-R.-pomeroyi-DSS-3-+-Pathways.tsv", sep="\t"):
    PATHWAYS = pd.read_csv(pathways_data, sep="\t")
    PATHWAYS_DICT = {r: p for i, (r, _, p) in PATHWAYS[PATHWAYS["Names"].notnull()].iterrows()}

    pathways = PATHWAYS_DICT.get(reaction_id, "").split(" // ")
    return pathways if len(pathways) > 1 else pathways[0]
