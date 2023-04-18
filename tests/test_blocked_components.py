from cobra.io import read_sbml_model
from cobra.flux_analysis.reaction import assess
from cobra.flux_analysis import find_blocked_reactions


def test_blocked_biomass():
    model = read_sbml_model("clean_models/Rpom_05.xml")
    
    # Returns True if biomass can be produced, otherwise a
    # dictionary of {precursor: Status, products: Status}
    assert assess(model, "RPOM_provisional_biomass")


def main():
    test_blocked_biomass()


if __name__ == "__main__":
    main()
