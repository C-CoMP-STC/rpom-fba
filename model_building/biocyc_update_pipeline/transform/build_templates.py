import os
import pickle
import pandas as pd
import warnings
from cobra.io import read_sbml_model
from tqdm import tqdm

# Utility function for clean retrieval of slot values


def wrap(slot_value):
    if not isinstance(slot_value, list):
        slot_value = [slot_value]
    return slot_value


def get_gene_reaction_rule(reaction_data, proteins):
    # Get enzymatic reactions associated with the reaction
    enzymatic_reactions = wrap(reaction_data.get(
        "enzymatic-reaction", {}).get("Enzymatic-Reaction", []))

    # Build rule, taking the "or" of all enzymatic reactions
    rule = []
    for enzymatic_reaction in enzymatic_reactions:
        # enzymatic_reaction = enzymatic_reaction["Enzymatic-Reaction"]

        # Gene-reaction-rule for this enzymatic reaction
        # (elements will be "and"ed together, and all sub-rules will be "or"ed together at the end)
        sub_rule = []

        # Get enzymes of enzymatic reaction
        enzymes = wrap(enzymatic_reaction.get("enzyme", []))
        for enzyme in enzymes:
            # Retrieve protein data
            enzyme = enzyme["Protein"]
            enzyme = proteins[enzyme["@frameid"]]

            # Check if the enzyme is a protein complex, by checking if it has components
            is_complex = "component" in enzyme

            # Get protein components of complexes
            if is_complex:
                components = [component["Protein"]
                              for component in wrap(enzyme["component"])]
            else:
                components = [enzyme]

            # Get gene frameids for each component
            for component in components:
                genes = [gene["Gene"]
                         for gene in wrap(component.get("gene", []))]
                gene_ids = [gene["@frameid"] for gene in genes]
                assert len(genes) <= 1, f"Expected <= 1 gene per polypeptide, got {len(genes)}"

                # Extend the sub-rule with the gene frameid
                sub_rule.extend(gene_ids)

        # Add the sub-rule to the rule, if it's not empty
        if sub_rule:
            rule.append(sub_rule)

    # Build rule as string
    rule_str = [(f"({' and '.join(sub_rule)})"
                 if len(sub_rule) > 1
                 else sub_rule[0])
                for sub_rule in rule]
    rule_str = f"({' or '.join(rule_str)})"

    # Return computed gene-reaction rule
    return rule_str


def get_stoichiometry(reaction_data):
    # Utility function for retrieving stoichiometry from reaction data dictionary
    # (as obtained from biocyc API request)
    stoich = {}

    # Get reactants and products
    left = wrap(reaction_data.get("left", None))
    right = wrap(reaction_data.get("right", None))

    for side, sign in [(left, -1), (right, 1)]:
        for compound in side:
            if compound is None:
                continue

            # Get data for the compound, which may be one of several forms (based on my observations).
            # This is so unnecessarily complicated! I've seen the following forms:
            #
            # 1. A bare string, indicating a compound name, e.g. "a supercoiled duplex DNA"
            #    (probably not matching to any actual object in the database)
            #
            # 2. A dictionary with some of the following keys in some combination (together with example values):
            #        - "Compound"/"Protein"/"RNA" : {"@frameid": ..., ...}
            #               The most "standard" key, representing the schema class of the compound,
            #               and containing the frameid and other standard data.
            #
            #        - "#text" : "a deoxyribonucleic acid"
            #              Typically used to name compounds that are not actually objects in the database,
            #              e.g. electrons in redox half-reactions ("e<SUP>-</SUP>") or generic compounds ("a deoxyribonucleic acid", "an antibiotic", ...).
            #              Use this as a fallback for ID if the frameid is not present.
            #
            #        - "name-slot" : "N-Name" / "N-1-Name" / "N+1-Name" / ...,
            #              Used in polymerization reactions to indicate the length of the polymer.
            #              I noticed at least one instance of "N-2-NAME", with the "NAME" capitalized, so take note and
            #              consider case-insensitivity. I'll tack this on to the end of the compound id, whether from @frameid or #text.
            #
            #        - "coefficient" : {"@datatype": "integer", "#text": "3"} or {"@datatype": "string", "#text": "n-1"}
            #               The stoichiometric coefficient of the compound (always positive, so we need to add our own signs).
            #               May also be a string like "n" or "n-1", representing an undefined number of compounds
            #               (presumably in reaction classes for polymerization). We'll preserve such values, but recommend deletion down the line.
            #               If not present, defaults to 1.
            #
            #        - "compartment" : {"cco" : {..., "@frameid": ...}},
            #              The compartment of the compound, as a dictionary with the frameid of the compartment.
            #              Many of these will be the abstract compartments "CCO-IN" and "CCO-OUT", which really should be
            #              clarified through a "RXN-LOCATIONS" slot in the PGDB, though that slot does not appear to be accessible like this.
            #              We'll assume that CCO-IN is the cytoplasm and CCO-OUT is the periplasm (represented in our model as [c] and [p], respectively).
            #              The default is the cytoplasm ('CCO-CYTOSOL').

            compound_id = ""
            coefficient = sign  # Default coefficient is -1 for reactants, 1 for products
            if isinstance(compound, str):
                # Reactant is a compound name, so we'll use the name as the ID
                # and assume it's in the cytoplasm
                compound_id = f"{compound}[c]"
            else:
                # Reactant is a dictionary, so we'll try to extract the frameid and coefficient
                # First, look for a compound, protein, or RNA:
                found_id = False
                for key in ["Compound", "Protein", "RNA"]:
                    if key in compound:
                        compound_id = compound[key]["@frameid"]
                        found_id = True
                        break

                if not found_id:
                    # If we didn't find a frameid, use the #text field
                    compound_id = compound.get("#text", "")

                # For polymerization reactions, add the polymer length to the end of the compound id
                if "name-slot" in compound:
                    compound_id += f"__{compound['name-slot']}"

                # Get the compartment
                compartment = compound.get("compartment", {}).get(
                    "cco", {}).get("@frameid", "CCO-CYTOSOL")
                if compartment == "CCO-CYTOSOL":
                    compartment = "c"
                elif compartment == "CCO-IN":
                    compartment = "c"
                elif compartment == "CCO-OUT":
                    compartment = "p"
                else:
                    warnings.warn(
                        f"Unknown compartment frameid: {compartment}")

                compound_id += f"[{compartment}]"

                # Get the stoichiometric coefficient (defaults to -1, see above)
                # and make it negative to indicate a reactant
                if "coefficient" in compound:
                    if compound["coefficient"]["@datatype"] == "integer":
                        coefficient = sign * \
                            int(compound["coefficient"]["#text"])
                    else:
                        # coefficient is a string, like "n" or "n-1"
                        coefficient = f"{'-' if sign == -1 else '+'}{compound['coefficient']['#text']}"

            # Store coefficient
            if isinstance(coefficient, int):
                # Sometimes, compounds appear on both sides,
                # so we need to take the net
                stoich[compound_id] = stoich.get(compound_id, 0) + coefficient
            else:
                # coefficient is a string
                stoich[compound_id] = stoich.get(compound_id, "") + coefficient

    # Final cleanup: remove any zero-coefficient compounds
    stoich = {met: coeff for met, coeff in stoich.items() if coeff != 0}

    return stoich


def get_bounds(reaction_data):
    bounds = None
    match reaction_data.get("reaction-direction", None):
        case "REVERSIBLE":
            bounds = (-1000, 1000)
        case "PHYSIOL-LEFT-TO-RIGHT" | "IRREVERSIBLE-LEFT-TO-RIGHT" | "LEFT-TO-RIGHT":
            bounds = (0, 1000)
        case "PHYSIOL-RIGHT-TO-LEFT" | "IRREVERSIBLE-RIGHT-TO-LEFT" | "RIGHT-TO-LEFT":
            bounds = (-1000, 0)
        case None:
            bounds = None
            # TODO: Compute based on enzymatic reactions as described at https://www.biocyc.org/PGDBConceptsGuide.shtml#TAG:__tex2page_sec_4.2
            warnings.warn(
                "No reaction direction specified (inference from enzymatic reactions not yet implemented)")
        case _:
            warnings.warn("Unknown reaction direction: " +
                          f"{reaction_data.get('reaction-direction', None)}")
    return bounds


def build_existing_reactions_template(model,
                                      reactions1,
                                      reactions2,
                                      proteins1,
                                      proteins2,
                                      substrate_exchanges=[
                                          "EX_glc",
                                          "EX_ac"
                                      ],
                                      skip_reactions=["ATPM", "Rpom_hwa_biomass"]):

    # Get flux vectors for the specified substrates,
    # to ensure that we keep reactions necessary for growth on these substrates
    flux_vectors = {}
    for ex in substrate_exchanges:
        with model:
            exchange = model.reactions.get_by_id(ex)
            exchange.bounds = (-10, 0)
            sol = model.optimize()

            if sol.objective_value > 0:
                flux_vectors[ex] = sol.fluxes

    report = []
    for reaction in tqdm(model.reactions, "Processing existing reactions..."):
        # Skip specified reactions, as well as any exchanges
        if reaction.id in skip_reactions or reaction in model.exchanges:
            continue

        line = {}

        # Store reaction ID
        line["Reaction ID"] = reaction.id

        # Store whether the reaction is an exchange/transport reaction, as we'll keep all of these
        line["Is-Exchange"] = reaction in model.exchanges
        line["Is-Transport"] = reaction.id.endswith("tex") or reaction.id.endswith("tpp")

        # Test whether the reaction has nonzero flux on specified substratess
        used_on_any = False
        for ex in flux_vectors.keys():
            used_on_ex = (flux_vectors[ex][reaction.id] != 0)
            line[f"Used on {ex}?"] = used_on_ex
            used_on_any |= used_on_ex
        line["Used (any)?"] = used_on_any

        # Check if the reaction is present in either biocyc database,
        # and retrieve data if so.
        reaction_db1 = reactions1.get(reaction.notes.get("stem", reaction.id), None)
        reaction_db2 = reactions2.get(reaction.notes.get("stem", reaction.id), None)
        line["In DB1?"] = reaction_db1 is not None
        line["In DB2?"] = reaction_db2 is not None

        # Compute and store gene-reaction-rules
        line["Gene-reaction rule 1"] = (get_gene_reaction_rule(reaction_db1, proteins1)
                                        if reaction_db1 is not None
                                        else None)
        line["Gene-reaction rule 2"] = (get_gene_reaction_rule(reaction_db2, proteins2)
                                        if reaction_db2 is not None
                                        else None)

        # Compile potential issues
        issues = []
        if not used_on_any:
            issues.append("Reaction is never used on the given substrates")
        if reaction_db1 is None and reaction_db2 is None:
            issues.append("Reaction is found in neither database")
        line["Issues"] = issues

        # Compute recommendation
        recommendation = None
        if (not used_on_any
            and not line["Is-Exchange"]
            and not line["Is-Transport"]
            and reaction_db1 is None
            and reaction_db2 is None):
            recommendation = "Delete"
        else:
            recommendation = "Keep"
        line["Recommendation"] = recommendation

        # Add line to report
        report.append(line)

    # Turn into a DataFrame and return
    report = pd.DataFrame(report)
    return report


def build_new_reactions_template(
        model,
        reactions1,
        reactions2,
        metabolites1,
        metabolites2,
        proteins1,
        proteins2):
    report = []

    # Look only at reactions not in model
    reactions_db1 = set(reactions1.keys())
    reactions_db2 = set(reactions2.keys())
    reactions_model = set([rxn.notes.get("stem", rxn.id) for rxn in model.reactions])
    reactions_to_add = (reactions_db1 | reactions_db2) - reactions_model

    for reaction in reactions_to_add:
        line = {}

        # Store reaction ID
        line["Reaction ID"] = reaction

        # Check if reaction is in either database, and pull the data
        in_db1 = reaction in reactions_db1
        in_db2 = reaction in reactions_db2
        reaction_data_1 = reactions1.get(reaction, None)
        reaction_data_2 = reactions2.get(reaction, None)

        # For some fields, favor data from DB2 if available
        # reaction_data = reaction_data_2 if in_db2 else reaction_data_1

        # For some fields, favor data from DB1 if available
        reaction_data = reaction_data_1 if in_db1 else reaction_data_2

        # Store reaction name, and the databases in which it is present
        line["Reaction name 1"] = (reaction_data_1.get(
            "common-name", {}).get("#text", None)
            if reaction_data_1 is not None
            else None)
        line["Reaction name 2"] = (reaction_data_2.get(
            "common-name", {}).get("#text", None)
            if reaction_data_2 is not None
            else None)
        line["In DB1?"] = in_db1
        line["In DB2?"] = in_db2

        # Compute gene-reaction rules
        line["Gene-reaction rule 1"] = (get_gene_reaction_rule(
            reaction_data_1, proteins1)
            if reaction_data_1 is not None
            else None)
        line["Gene-reaction rule 2"] = (get_gene_reaction_rule(
            reaction_data_2, proteins2)
            if reaction_data_2 is not None
            else None)

        # Stoichiometry
        stoichiometry_1 = get_stoichiometry(
            reaction_data_1) if in_db1 else None
        stoichiometry_2 = get_stoichiometry(
            reaction_data_2) if in_db2 else None
        line["Stoichiometry 1"] = stoichiometry_1
        line["Stoichiometry 2"] = stoichiometry_2

        # Prefer stoichiometry from DB2, if available
        # stoichiometry = stoichiometry_2 if in_db2 else stoichiometry_1

        # Prefer stoichiometry from DB1, if available
        stoichiometry = stoichiometry_1 if in_db1 else stoichiometry_2

        # Get bounds
        bounds_1 = get_bounds(reaction_data_1) if in_db1 else None
        bounds_2 = get_bounds(reaction_data_2) if in_db2 else None
        line["Bounds 1"] = bounds_1
        line["Bounds 2"] = bounds_2

        # Check if any metabolites in the stoichiometry are missing from the model
        missing_metabolites = [
            met for met in stoichiometry.keys() if met not in model.metabolites]
        line["Metabolites not in Model"] = missing_metabolites

        # Check if any missing metabolites are also missing from the databases
        # (this can happen in the case of "pseudo-metabolites" like "a deoxyribonucleic acid")
        fake_metabolites = [met
                            for met in missing_metabolites
                            if met.split('[')[0] not in metabolites1
                            and met.split('[')[0] not in metabolites2]
        line["Metabolites without data"] = fake_metabolites

        # Check for class metabolites
        # metabolites = metabolites2 if in_db2 else metabolites1
        metabolites = metabolites1 if in_db1 else metabolites2
        has_class_metabolites = False
        for met in stoichiometry.keys():
            met = met.split("[")[0]
            if metabolites.get(met, {}).get("@class", False) == "true":
                has_class_metabolites = True
                break
        line["Has-Class-Metabolites"] = has_class_metabolites

        # Check metabolites for generic polymers (...__N<>)
        polymerization = any("__N" in met for met in stoichiometry.keys())
        line["Polymerization"] = polymerization

        # Check if reaction is class
        line["Is-Class"] = reaction_data.get("@class", False) == "true"

        # Check if reaction has subclasses or instances
        line["Has-Subclasses"] = reaction_data.get(
            "subclass", None) is not None
        line["Has-Instances"] = reaction_data.get("instance", None) is not None

        # Check for problematic compartments ([p], [CCO-MIDDLE], etc.)
        compartments_to_review = {met[met.find(
            "["):] for met in stoichiometry.keys() if met[met.find("["):] != "[c]"}

        # Compile potential issues
        issues = []
        has_fake_metabolites = (len(fake_metabolites) > 0)
        no_stoich = (stoichiometry is None or len(stoichiometry) == 0)
        stoich_disagreement = (
            stoichiometry_1 is not None and stoichiometry_2 is not None and stoichiometry_1 != stoichiometry_2)
        class_with_instances = line["Is-Class"] and line["Has-Instances"]
        bounds_disagree = (
            bounds_1 is not None and bounds_2 is not None and bounds_1 != bounds_2)

        if no_stoich:
            issues.append("No stoichiometry data")
        if stoich_disagreement:
            issues.append("Stoichiometry differs across databases")
        if missing_metabolites:
            issues.append("Reaction has metabolites missing from the model")
        if has_fake_metabolites:
            issues.append("Reaction has metabolites without data")
        if has_class_metabolites:
            issues.append(
                "Reaction has class metabolites, may require specialization")
        if polymerization:
            issues.append("Reaction is a polymerization reaction")
        if class_with_instances:
            issues.append("Reaction is a class reaction with instances")
        if bounds_disagree:
            issues.append("Bounds disagree across databases")
        if compartments_to_review:
            issues.append(
                "Reaction has metabolites in non-cytoplasm compartments")
        line["Issues"] = issues

        # Compute recommendation
        recommendation = ""
        if no_stoich or class_with_instances or has_fake_metabolites:
            recommendation = "Do not add"
        elif stoich_disagreement or has_class_metabolites or polymerization or bounds_disagree or compartments_to_review:
            recommendation = "Manual review"
        else:
            # May have "missing" metabolites, but those metabolites
            # can be added into the model
            recommendation = "Add"
        line["Recommendation"] = recommendation

        # Add line to report
        report.append(line)

    return pd.DataFrame(report)


def build_new_metabolites_template(model, candidate_mets, metabolites1, metabolites2):
    report = []
    for met in candidate_mets:
        line = {}

        # Check if metabolite is in either database
        in_db1 = met in metabolites1
        in_db2 = met in metabolites2
        if not in_db1 and not in_db2:
            continue

        # Prefer data from DB2, if available
        # met_data = metabolites2[met] if in_db2 else metabolites1[met]

        # Prefer data from DB1, if available
        met_data = metabolites1[met] if in_db1 else metabolites2[met]

        # Store metabolite ID
        line["Metabolite ID"] = met
        line["Metabolite Name"] = met_data.get(
            "common-name", {}).get("#text", None)

        # Check if metabolite is in the model already
        line["In Model?"] = (f"{met}[c]" in model.metabolites
                             or f"{met}[p]" in model.metabolites
                             or f"{met}[e]" in model.metabolites)

        # Get formula, charge
        formula = met_data.get("cml", {}).get("molecule", {}).get(
            "formula", {}).get("@concise", None)
        line["Formula"] = "".join(formula.split()) if isinstance(formula, str) else formula

        charge = met_data.get("cml", {}).get(
            "molecule", {}).get("@formalCharge", None)
        charge = int(charge) if charge is not None else None
        line["Charge"] = charge

        # Check if metabolite is a class
        line["Is-Class"] = met_data.get("@class", False) == "true"

        # Check if metabolite has subclasses or instances
        instances = [c for c in met_data["instance"]
                     ] if "instance" in met else None
        subclasses = [c for c in met_data["subclass"]
                      ] if "subclass" in met else None
        line["Instances"] = instances
        line["Subclasses"] = subclasses

        # Add line to report
        report.append(line)

    return pd.DataFrame(report)


def build_new_genes_template(genes1, genes2):
    report = []

    # Get all gene frameids
    gene_frameids = set(genes1.keys()) | set(genes2.keys())
    for gene in gene_frameids:
        line = {}

        # Check if gene is in either database
        in_db1 = gene in genes1
        in_db2 = gene in genes2

        # Prefer data from DB2, if available
        # gene_data = genes2[gene] if in_db2 else genes1[gene]

        # Prefer data from DB1, if available
        gene_data = genes1[gene] if in_db1 else genes2[gene]

        # Store gene ID, name, synonyms
        line["Gene ID"] = gene
        line["Gene Name"] = gene_data.get("common-name", {}).get("#text", None)
        synonyms = ([syn["#text"] for syn in wrap(gene_data["synonym"])]
                    if "synonym" in gene_data
                    else None)
        line["Synonyms"] = synonyms

        # Store database(s) of origin
        line["In DB1?"] = in_db1
        line["In DB2?"] = in_db2

        # Store position
        left = gene_data.get("left-end-position", {}).get("#text", None)
        right = gene_data.get("right-end-position", {}).get("#text", None)
        line["Left-Position"] = int(left) if left is not None else None
        line["Right-Position"] = int(right) if right is not None else None

        # Store replicon
        replicon = gene_data["replicon"]["Genetic-Element"]["@frameid"]
        line["Replicon"] = replicon

        report.append(line)

    return pd.DataFrame(report)


def main():
    DB1 = "RUEGERIA_POMEROYI_DSS3"
    DB2 = "GCF_000011965"
    MODEL = "model/Rpom_05.xml"

    # Load current model
    model = read_sbml_model(MODEL)

    # Load raw data from files =========================================================================
    DATA_DIR = "model_building/biocyc_update_pipeline/data/"

    # Reactions
    with open(os.path.join(DATA_DIR, f"raw_reaction_data__{DB1}.pkl"), "rb") as f:
        reactions1 = pickle.load(f)
    with open(os.path.join(DATA_DIR, f"raw_reaction_data__{DB2}.pkl"), "rb") as f:
        reactions2 = pickle.load(f)
    # Reformat reaction data to be keyed by id
    reactions1 = {rxn["@frameid"]: rxn for rxn in reactions1}
    reactions2 = {rxn["@frameid"]: rxn for rxn in reactions2}

    # Metabolites
    with open(os.path.join(DATA_DIR, f"raw_metabolite_data__{DB1}.pkl"), "rb") as f:
        metabolites1 = pickle.load(f)
    with open(os.path.join(DATA_DIR, f"raw_metabolite_data__{DB2}.pkl"), "rb") as f:
        metabolites2 = pickle.load(f)
    # Reformat metabolite data to be keyed by id
    # TODO: Using a hack to compensate for failure to unwrap RNA metabolites
    # at the extract stage. Remove this once the fixed extract stage is run.
    metabolites1 = {
        (
            met["@frameid"]
            if "@frameid" in met
            else met["RNA"]["@frameid"]
        ): (met if "@frameid" in met else met["RNA"])
        for met in metabolites1}
    metabolites2 = {
        (
            met["@frameid"]
            if "@frameid" in met
            else met["RNA"]["@frameid"]
        ): (met if "@frameid" in met else met["RNA"])
        for met in metabolites2}

    # Proteins
    with open(os.path.join(DATA_DIR, f"raw_protein_data__{DB1}.pkl"), "rb") as f:
        proteins1 = pickle.load(f)
    with open(os.path.join(DATA_DIR, f"raw_protein_data__{DB2}.pkl"), "rb") as f:
        proteins2 = pickle.load(f)
    proteins1 = {protein["@frameid"]: protein for protein in proteins1}
    proteins2 = {protein["@frameid"]: protein for protein in proteins2}

    # Genes
    with open(os.path.join(DATA_DIR, f"raw_gene_data__{DB1}.pkl"), "rb") as f:
        genes1 = pickle.load(f)
    with open(os.path.join(DATA_DIR, f"raw_gene_data__{DB2}.pkl"), "rb") as f:
        genes2 = pickle.load(f)
    genes1 = {gene["@frameid"]: gene for gene in genes1}
    genes2 = {gene["@frameid"]: gene for gene in genes2}

    # Build templates ==================================================================================

    # Reactions already in the model -
    # update gene-reaction rules, check for usage on key substrates,
    # and recommend deletion if not used
    existing_reactions_template = build_existing_reactions_template(
        model, reactions1, reactions2, proteins1, proteins2,
        substrate_exchanges=["EX_glc",
                             "EX_ac",
                             "EX_CPD0-1265",
                             "EX_CPD-335",
                             "EX_CARNITINE",
                             "EX_CHOLINE",
                             "EX_CIT",
                             "EX_L-CYSTEATE",
                             "EX_CPD-12693",
                             "EX_SS-DIMETHYL-BETA-PROPIOTHETIN",
                             "EX_ECTOINE",
                             "EX_FUM",
                             "EX_MAL",
                             "EX_SUC",
                             "EX_BETA-D-XYLOSE",
                             "EX_GLYCEROL",
                             "EX_GLYCEROL-3P",
                             "EX_CPD-3745",
                             "EX_N-ACETYL-D-GLUCOSAMINE",
                             "EX_SPERMIDINE",
                             "EX_CADAVERINE",
                             "EX_PUTRESCINE",
                             "EX_TAURINE",
                             "EX_THYMIDINE",
                             "EX_TRIMETHYLAMINE-N-O",
                             "EX_VAL",
                             "EX_GLT"])

    # Candidate new reactions to add -
    # check for missing metabolites, class metabolites, polymerization reactions,
    # and recommend manual review for these cases
    new_reactions_template = build_new_reactions_template(
        model, reactions1, reactions2, metabolites1, metabolites2, proteins1, proteins2)

    # Get all candidate new metabolites to add
    all_new_mets = set()
    for _, (stoich1, stoich2) in new_reactions_template[["Stoichiometry 1", "Stoichiometry 2"]].iterrows():
        all_new_mets.update(stoich1.keys() if stoich1 is not None else set())
        all_new_mets.update(stoich2.keys() if stoich2 is not None else set())
    # Remove compartments
    all_new_mets = {met.split("[")[0] for met in all_new_mets}

    # Build new metabolites template
    new_metabolites_template = build_new_metabolites_template(
        model, all_new_mets, metabolites1, metabolites2)

    # Build new genes template
    new_genes_template = build_new_genes_template(genes1, genes2)

    # Ensure file exists
    OUT = "model_building/biocyc_update_pipeline/templates/template.xlsx"
    if not os.path.exists(OUT):
        pd.DataFrame().to_excel(OUT)

    # Write in overwrite mode
    with pd.ExcelWriter(
            OUT,
            mode="w") as writer:
        existing_reactions_template.to_excel(writer,
                                             index=False,
                                             sheet_name="Existing reactions")
        new_reactions_template.to_excel(writer,
                                        index=False,
                                        sheet_name="New reactions")
        new_metabolites_template.to_excel(writer,
                                          index=False,
                                          sheet_name="New metabolites")
        new_genes_template.to_excel(writer,
                                    index=False,
                                    sheet_name="New genes")


if __name__ == "__main__":
    main()
