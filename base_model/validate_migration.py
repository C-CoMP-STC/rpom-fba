from cobra.io import validate_sbml_model

def validate_model_migration(sbml_model, sbml_model_file, reversibility):
    validate_reversibility(sbml_model, reversibility)
    _, errors = validate_sbml_model(sbml_model_file)
    for error_type, error_list in errors.items():
        print(f"{error_type}: {len(error_list)} errors")

        if len(error_list) <= 10:
            for err in error_list:
                print(f"\t{err}")
    

def validate_reversibility(sbml_model, reversibility):
    for rxn, rev in zip(sbml_model.reactions, reversibility):
        lb, _ = rxn.bounds
        if rev and lb >= 0:
            print(f"Reaction {rxn} expected to be reversible, but lower bound is {lb} >= 0.")
        elif not rev and lb < 0:
            print(f"Reaction {rxn} expected to be irreversible, but lower bound is {lb} < 0.")