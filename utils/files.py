import os
import pandas as pd
from cobra.io import read_sbml_model

def model_to_excel(model, outfile):
    metabolites = {
        "id": [],
        "name": [],
        "annotation" : [],
        "charge" : [],
        "compartment" : [],
        "elements": [],
        "formula": [],
        "formula_weight": [],
        "notes": [],
        "reactions": []
    }
    for met in model.metabolites:
        metabolites["id"].append(met.id)
        metabolites["name"].append(met.name)
        metabolites["annotation"].append(met.annotation)
        metabolites["charge"].append(met.charge)
        metabolites["compartment"].append(met.compartment)
        metabolites["elements"].append(met.elements)
        metabolites["formula"].append(met.formula)
        metabolites["formula_weight"].append(met.formula_weight)
        metabolites["notes"].append(met.notes)
        metabolites["reactions"].append("{" + ", ".join(rxn.id for rxn in met.reactions) + "}")
    metabolites = pd.DataFrame(metabolites)

    genes = {
        "id" : [],
        "name": [],
        "annotation": [],
        "functional": [],
        "notes": [],
        "reactions": [],
    }
    for gene in model.genes:
        genes["id"].append(gene.id)
        genes["name"].append(gene.name)
        genes["annotation"].append(gene.annotation)
        genes["functional"].append(gene.functional)
        genes["notes"].append(gene.notes)
        genes["reactions"].append("{" + ", ".join(rxn.id for rxn in gene.reactions) + "}")
    genes = pd.DataFrame(genes)

    reactions = {
        "id": [],
        "name": [],
        "annotation": [],
        "boundary": [],
        "lower_bound": [],
        "upper_bound": [],
        "compartments": [],
        "functional": [],
        "gene_reaction_rule": [],
        "genes": [],
        "metabolites": [],
        "products": [],
        "reactants": [],
        "reaction": [],
        "reverse_id": [],
        "reversibility": []
    }
    for rxn in model.reactions:
        reactions["id"].append(rxn.id)
        reactions["name"].append(rxn.name)
        reactions["annotation"].append(rxn.annotation)
        reactions["boundary"].append(rxn.boundary)
        reactions["lower_bound"].append(rxn.lower_bound)
        reactions["upper_bound"].append(rxn.upper_bound)
        reactions["compartments"].append(rxn.compartments)
        reactions["functional"].append(rxn.functional)
        reactions["gene_reaction_rule"].append(rxn.gene_reaction_rule)
        reactions["genes"].append("{" +
                                  ", ".join(gene.id for gene in rxn.genes) +
                                  "}")
        reactions["metabolites"].append("{" +
                                        ", ".join(f"{met} : {coeff}" for met, coeff in rxn.metabolites.items()) +
                                        "}")
        reactions["products"].append([met.id for met in rxn.products])
        reactions["reactants"].append([met.id for met in rxn.reactants])
        reactions["reaction"].append(rxn.reaction)
        reactions["reverse_id"].append(rxn.reverse_id)
        reactions["reversibility"].append(rxn.reversibility)
    reactions = pd.DataFrame(reactions)


    with pd.ExcelWriter(outfile) as writer:
        metabolites.to_excel(writer, sheet_name='Metabolites', index=False)
        genes.to_excel(writer, sheet_name="Genes", index=False)
        reactions.to_excel(writer, sheet_name="Reactions", index=False)


def main():
    MODEL = "model/Rpom_05.xml"
    OUTFILE = "out/model/Rpom_05.xlsx"

    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

    model = read_sbml_model(MODEL)
    model_to_excel(model, OUTFILE)

if __name__ == "__main__":
    main()
