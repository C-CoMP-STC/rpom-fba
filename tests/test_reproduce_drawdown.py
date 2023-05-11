import pytest
import json
from cobra.io import read_sbml_model

from tests.testing_utils import TEST_MODEL


def get_carbon_sources(carbon_source_ids="parameters/uptake_rates/carbon_source_ids.json"):
    with open(carbon_source_ids, "r") as f:
        carbon_source_ids = json.load(f)

    return list(carbon_source_ids.values())


@pytest.mark.parametrize("metabolite", get_carbon_sources())
def test_growth_and_uptake(metabolite):
    model = read_sbml_model(TEST_MODEL)


def main():
    for carbon_source in get_carbon_sources():
        test_growth_and_uptake(carbon_source)


if __name__ == "__main__":
    main()
