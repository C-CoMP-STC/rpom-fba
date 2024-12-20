import pytest
from cobra import Reaction
from cobra.flux_analysis.reaction import assess
from cobra.io import read_sbml_model

from tests.testing_utils import TEST_MODEL, run_all_tests_in_object
from utils.cobra_utils import set_active_bound


class TestBlockedComponents():
    model = read_sbml_model(TEST_MODEL)

    ex_glc = model.reactions.get_by_id("EX_glc")
    set_active_bound(ex_glc, abs(float(ex_glc.annotation["Experimental rate"])))

    def test_blocked_biomass(self):
        # Returns True if biomass can be produced, otherwise a
        # dictionary of {precursor: Status, products: Status}
        can_produce_biomass = assess(self.model, "RPOM_provisional_biomass")
        message = ""
        if isinstance(can_produce_biomass, dict):
            if isinstance(can_produce_biomass['precursors'], dict):
                blocked_precursors = '\n'.join(
                    met.id for met in can_produce_biomass['precursors'])
                message += f"Biomass blocked by production of:\n {blocked_precursors}"
            if isinstance(can_produce_biomass['products'], dict):
                non_absorbed_products = '\n'.join(
                    met.id for met in can_produce_biomass['products'])
                message += f"\n Biomass blocked by absorption of:\n {non_absorbed_products}"
        assert can_produce_biomass == True, message

    @pytest.mark.parametrize("component,coefficient",
                             (model
                              .reactions
                              .get_by_id("RPOM_provisional_biomass")
                              .metabolites
                              .items()))
    def test_biomass_component(self, component, coefficient):
        with self.model as model:
            model.reactions.get_by_id(
                "RPOM_provisional_biomass").bounds = (0, 0)

            single_metabolite_objective = Reaction(f"{component.id}_objective")
            single_metabolite_objective.add_metabolites({
                component: coefficient
            })
            model.add_reactions([single_metabolite_objective])

            model.objective = single_metabolite_objective
            solution = model.optimize()

            assert solution.objective_value > 0, (f"Unable to {'produce' if coefficient < 0 else 'absorb'} "
                                                f"{component} (objective value = {solution.objective_value} < {coefficient})")


def main():
    tests = TestBlockedComponents()
    run_all_tests_in_object(tests,
                            test_params={
                                "test_biomass_component": tests.model.reactions.get_by_id("RPOM_provisional_biomass").metabolites.items()
                            })


if __name__ == "__main__":
    main()
