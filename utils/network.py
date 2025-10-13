import networkx as nx

def model_to_network(model):
    g = nx.Graph()

    # Add nodes for metabolites, reactions
    g.add_nodes_from([
        (met.id,
         {
             "node_type": "Metabolite",
             "name": met.name,
             "formula" : met.formula,
             "charge" : met.charge,
             "compartment" : met.compartment
        })
        for met in model.metabolites
    ])

    g.add_nodes_from([
        (rxn.id,
         {
             "node_type": "Reaction",
             "name": rxn.name,
             "reaction" : rxn.reaction,
             "gpr" : rxn.gene_reaction_rule,
             "lower" : rxn.lower_bound,
             "lower" : rxn.upper_bound,
             "boundary" : rxn in model.boundary
        })
        for rxn in model.reactions
    ])

    # Add edges between reactions and metabolites
    for rxn in model.reactions:
        g.add_edges_from([
            (met.id, rxn.id,
             {"coefficient": coeff})
            for met, coeff in rxn.metabolites.items()
        ])
    
    return g


if __name__ == "__main__":
    from cobra.io import read_sbml_model

    model = read_sbml_model("model/Rpom_05.xml")
    g = model_to_network(model)