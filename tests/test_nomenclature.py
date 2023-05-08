from cobra.io import read_sbml_model
from cobra.manipulation.validate import check_metabolite_compartment_formula

from tests.testing_utils import TEST_MODEL, run_all_tests_in_object


class TestNomenclature:
    model = read_sbml_model(TEST_MODEL)

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
