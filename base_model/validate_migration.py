from cobra.io import validate_sbml_model

def validate_model_migration(sbml_model_file):
    _, errors = validate_sbml_model(sbml_model_file)
    for error_type, error_list in errors.items():
        print(f"{error_type}: {len(error_list)} errors")

        if len(error_list) <= 10:
            for err in error_list:
                print(f"\t{err}")
