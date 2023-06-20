import json

import matplotlib
import pandas as pd
import pytest

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model

from tests.testing_utils import TEST_MODEL


CARBON_SOURCE_IDS = "parameters/uptake_rates/carbon_source_ids.json"
DRAWDOWN_DATA = "data/drawdown_clean.csv"
GROWTH_DATA = "data/growth_curves_raw.csv"


def get_carbon_sources(carbon_source_ids=CARBON_SOURCE_IDS):
    with open(carbon_source_ids, "r") as f:
        carbon_source_ids = json.load(f)

    return list(carbon_source_ids.values())


@pytest.fixture
def drawdown_data(drawdown_data=DRAWDOWN_DATA):
    return pd.read_csv(DRAWDOWN_DATA)

@pytest.fixture
def testing_model():
    return read_sbml_model(TEST_MODEL)


@pytest.mark.parametrize("metabolite", get_carbon_sources())
def test_growth_and_uptake(metabolite, drawdown_data, testing_model):

    # Set media to minimal media
    model = testing_model
    model.medium
    # run dFBA on minimal (seawater?) + metabolite


    # compare to growth data, drawdown data (what cutoff?)


    # plot differences
    


def main():
    for carbon_source in get_carbon_sources():
        drawdown_data = pd.read_csv(DRAWDOWN_DATA)
        model = read_sbml_model(TEST_MODEL)
        test_growth_and_uptake(carbon_source, drawdown_data, model)


if __name__ == "__main__":
    main()
