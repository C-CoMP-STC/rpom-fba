import os
from collections import defaultdict

from tqdm import tqdm
import pandas as pd

from pathway_analysis.access_db import access_rpom_biocyc


PATHWAY_NAME_TEMPLATE = "PWY-*"


def main():
    GENES_OUT = "out/biocyc_model_comparison/genes.csv"

    rpom = access_rpom_biocyc()

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

            for gene_id in rpom.genes_of_reaction(reaction_id):
                # Get gene data
                gene = rpom[gene_id]

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

                result["Pathway"].append(pathway.common_name)
                result["Pathway synonyms"].append(", ".join(pathway.names))
                result["Pathway BioCyc Frame ID"].append(pathway.frameid)
                result["Gene name"].append(gene.common_name)
                result["Gene synonyms"].append(", ".join(gene.names))
                result["Gene BioCyc Frame ID"].append(gene.frameid)
                result["Reaction EC Number(s)"].append(", ".join(reaction.ec_number if reaction.ec_number is not None else []))
                result["Reversible in BioCyc"].append(REVERSIBILITY[reaction.reaction_direction])
                result["Reversible in model"].append(None)
                result["In Moran lab genome?"].append(None)
                result["In model?"].append(None)

    result = pd.DataFrame(result)

    os.makedirs(os.path.dirname(GENES_OUT), exist_ok=True)
    with open(GENES_OUT, "w") as f:
        result.to_csv(f, index=False)

    # Save additional metrics:
    # - Number of reactions not accounted for by pathways in BioCyc
    # - Reactions and genes unique to Biocyc
    # - Reactions and genes unique to model
    # - "Coverage" for each BioCyc reaction (percent of reactions in model)


if __name__ == "__main__":
    main()
