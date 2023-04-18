from cobra.io import read_sbml_model
from cobra.manipulation.validate import check_mass_balance, check_metabolite_compartment_formula


def main():
    model = read_sbml_model("clean_models/Rpom_05.xml")
    #TODO: make into a real test after removing exchanges and biomass
    check_mass_balance(model)

    assert len(check_metabolite_compartment_formula(model)) == 0


if __name__ == "__main__":
    main()
