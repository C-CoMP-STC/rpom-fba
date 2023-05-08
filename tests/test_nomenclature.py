from cobra.io import read_sbml_model
from cobra.manipulation.validate import check_metabolite_compartment_formula

from testing_utils import run_all_tests_in_object

class TestNomenclature:
    def __init__(self, model="clean_models/Rpom_05.xml"):
        self.model = read_sbml_model(model)

    def test_formulas(self):
        # Check all metabolite formulas are well-formed
        assert len(check_metabolite_compartment_formula(self.model)) == 0

    def test_boundary(self):
        # Check all demands, sinks, and exchanges are labeled with
        # appropriate prefixes.
        assert all(rxn.id[:2] == "DM" for rxn in self.model.demands)
        assert all(rxn.id[:2] == "SK" for rxn in self.model.sinks)
        assert all(rxn.id[:2] == "EX" for rxn in self.model.exchanges)


def main():
    tests = TestNomenclature()
    run_all_tests_in_object(tests)


if __name__ == "__main__":
    main()
