from cobra.io import read_sbml_model
from cobra.manipulation.validate import check_mass_balance

from testing_utils import run_all_tests_in_object

class TestMassBalance:
    def __init__(self, model="clean_models/Rpom_05.xml"):
        self.model = read_sbml_model(model)
    
    def test_mass_balance(self):
        # Get all unbalanced reactions, then exclude
        # 1) transport, and
        # 2) biomass reactions.
        # Any other unbalanced reactions are a problem!
        unbalanced_reactions = check_mass_balance(self.model)
        unbalanced_reactions = {k: v
                                for k, v in unbalanced_reactions.items()
                                if k not in self.model.boundary}
        unbalanced_reactions = {k: v
                                for k, v in unbalanced_reactions.items()
                                if "biomass" not in k.id.lower()}
        message = f"{','.join(rxn.id for rxn in unbalanced_reactions.keys())} not mass balanced!"
        assert len(unbalanced_reactions) == 0, message


def main():
    tests = TestMassBalance()
    run_all_tests_in_object(tests)


if __name__ == "__main__":
    main()
