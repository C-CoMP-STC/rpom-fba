from cobra.io import read_sbml_model
from cobra.flux_analysis.reaction import assess

from tests.testing_utils import TEST_MODEL


def test_blocked_biomass():
    model = read_sbml_model(TEST_MODEL)

    # Returns True if biomass can be produced, otherwise a
    # dictionary of {precursor: Status, products: Status}
    can_produce_biomass = assess(model, "RPOM_provisional_biomass")
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


def main():
    test_blocked_biomass()


if __name__ == "__main__":
    main()
