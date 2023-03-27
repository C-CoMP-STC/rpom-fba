import os
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from tqdm import tqdm

from pathway_analysis.access_db import access_rpom_biocyc

PATHWAY_NAME_TEMPLATE = "PWY-*"


def main():
    GENES_OUT = "out/biocyc_model_comparison/genes.csv"
    PATHWAYS_OUT = "out/biocyc_model_comparison/pathways.csv"
    ORDERED_PATHWAYS = "pathway_analysis/pathway_order.txt"

    # Get biocyc database and model
    rpom = access_rpom_biocyc()
    model = read_sbml_model("clean_models/Rpom_05.xml")

    # Create report comparing BioCyc with model
    result = defaultdict(list)

    for pathway in tqdm(rpom.pathways, desc="Comparing pathways...", total=len(rpom.pathways.instances)):
        for reaction_id in pathway.reaction_list:
            reaction = rpom[reaction_id].get_frame_data()

            # NOTE: Some reactions in pathways are not actual reactions, but rather,
            # references to other pathways that feed into/out of this one. These should be skipped.
            #
            # This check may not be the best way to do it - I'm not sure if there's a more general way to
            # get the BioCyc schema element of a PFrame - so replace it if you find a better way!
            if reaction.instance_name_template == PATHWAY_NAME_TEMPLATE:
                continue

            genes = rpom.genes_of_reaction(reaction_id)
            genes = genes if len(genes) >= 1 else [None]
            for gene_id in genes:

                pathway_synonyms = "; ".join(
                    set(pathway.names) - {pathway.common_name})

                # Get gene data if there is a gene
                gene_common_name = None
                gene_synonyms = None
                gene_frameid = None
                if gene_id:
                    gene = rpom[gene_id]
                    gene_common_name = gene.common_name
                    gene_synonyms = "; ".join(
                        set(gene.names) - {gene.common_name})
                    gene_frameid = gene.frameid

                # See the BioCyc schema documentation for explanation of these interpretations.
                REVERSIBILITY = {
                    "|REVERSIBLE|": "Yes",
                    "|LEFT-TO-RIGHT|": "Maybe",
                    "|RIGHT-TO-LEFT|": "Maybe",
                    "|PHYSIOL-LEFT-TO-RIGHT|": "Physiologically No",
                    "|PHYSIOL-RIGHT-TO-LEFT|": "Physiologically No",
                    "|IRREVERSIBLE-LEFT-TO-RIGHT|": "No",
                    "|IRREVERSIBLE-RIGHT-TO-LEFT|": "No"
                }

                # TODO: Quick and dirty removal of "bars" - fix this later
                reaction_name_friendly = reaction.frameid[1:-1]
                model_matches = model.reactions.query(reaction_name_friendly)
                in_model = "No"
                reversible_in_model = None
                if len(model_matches) == 1:
                    model_rxn = model_matches[0]
                    in_model = "Yes"
                    bounds_reversible = model_rxn.lower_bound < 0 and model_rxn.upper_bound > 0
                    labeled_reversible = model_rxn.reversibility

                    reversible_in_model = "Yes" if labeled_reversible else "No"
                    if bounds_reversible != labeled_reversible:
                        reversible_in_model = "Inconsistent!"

                elif len(model_matches) > 1:
                    in_model = "Maybe"

                result["Pathway"].append(pathway.common_name)
                result["Pathway synonyms"].append(pathway_synonyms)
                result["Pathway BioCyc Frame ID"].append(pathway.frameid)
                result["Reaction BioCyc ID"].append(reaction.frameid)
                result["Reaction LEFT"].append("; ".join(reaction.left))
                result["Reaction RIGHT"].append("; ".join(reaction.right))
                result["Reaction EC Number(s)"].append(
                    "; ".join(reaction.ec_number if reaction.ec_number is not None else []))
                result["Gene name"].append(gene_common_name)
                result["Gene synonyms"].append(gene_synonyms)
                result["Gene BioCyc Frame ID"].append(gene_frameid)
                result["Reversible in BioCyc"].append(
                    REVERSIBILITY[reaction.reaction_direction])
                result["Reversible in model"].append(reversible_in_model)
                result["In Moran lab genome?"].append(None)
                result["Reaction In model?"].append(in_model)

    result = pd.DataFrame(result)

    # If present, re-order output according to specified order
    try:
        with open(ORDERED_PATHWAYS, "r") as f:
            order = f.read().split("\n")

        result = result.sort_values(by="Pathway",
                                    ascending=True,
                                    kind="mergesort",
                                    key=lambda c: [
                                        order.index(p) if p in order else np.Infinity for p in c])
    except:
        pass

    os.makedirs(os.path.dirname(GENES_OUT), exist_ok=True)
    with open(GENES_OUT, "w") as f:
        # Google sheets gets confused and starts merging cells
        # if empty cells aren't filled with " "
        result.to_csv(f, index=False, na_rep=" ")

    # Save additional metrics:
    # - "Coverage" for each BioCyc reaction (percent of reactions in model)
    # - Number of reactions not accounted for by pathways in BioCyc
    # - Reactions and genes unique to Biocyc
    # - Reactions and genes unique to model
    with open(PATHWAYS_OUT, "w") as f:
        timestamp = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        (result.groupby("Pathway").aggregate(
            Coverage=pd.NamedAgg(column="Reaction In model?",
                                 aggfunc=lambda s: f"{100 * (s == 'Yes').sum() / len(s):.1f} %"),
            Reactions_in_Model=pd.NamedAgg(column="Reaction In model?",
                                           aggfunc=lambda s: (s == 'Yes').sum()),
            Reactions_in_BioCyc=pd.NamedAgg(column="Reaction In model?",
                                           aggfunc=lambda s: len(s)),
            Calculated=pd.NamedAgg(column="Reaction In model?",
                                   aggfunc=lambda s: timestamp))
            .reset_index()
            .sort_values(by="Coverage",
                         ascending=False,
                         kind="mergesort",
                         key=lambda s: [float(c.strip("%")) for c in s])
            .sort_values(by="Pathway",
                         ascending=True,
                         kind="mergesort",
                         key=lambda c: [
                             order.index(p) if p in order else np.Infinity for p in c])
         ).to_csv(f, index=False, na_rep=" ")


if __name__ == "__main__":
    main()
