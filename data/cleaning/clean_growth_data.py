import os
import json
import pickle
import numpy as np
import pandas as pd

from parameters.od_calibration import CellDensityRegressor, DRAWDOWN_SCALE_FACTOR
from parameters.drawdown import MASS_PER_CELL, COLONY_VOLUME
from utils.units import u
from data.files import GROWTH_DATA_RAW, CARBON_SOURCES, GROWTH_DATA_CLEAN


def main():
    data = pd.read_csv(GROWTH_DATA_RAW)

    with open("parameters/conversions/od_to_cell_density.pickle", "rb") as f:
        od_reg = pickle.load(f)

    with open(CARBON_SOURCES, "r") as f:
        carbon_sources = json.load(f).keys()

    # TODO: Subtract blank?
    blank = data["mean blank"].values[1:].astype(float)
    time = data["time (h)"].values[1:].astype(float)

    result = {}
    for carbon_source in carbon_sources:
        # Subset to data pertaining to this carbon source
        source_data = data.loc[:, [c.startswith(carbon_source) for c in data.columns]].values[1:].astype(float)

        # Take the mean across rows, multiplying by the necessary scale factor
        mean_growth = source_data.mean(axis=1) * DRAWDOWN_SCALE_FACTOR
        result[f"{carbon_source}_mean_OD"] = mean_growth

        # Translate to cell count
        cell_count = (od_reg.predict(mean_growth[~np.isnan(mean_growth)]).T[0] /u.mL) * COLONY_VOLUME
        result[f"{carbon_source}_predicted_count"] = cell_count.to("dimensionless").magnitude

        # Translate to mass
        total_mass = MASS_PER_CELL * cell_count
        result[f"{carbon_source}_predicted_mass"] = total_mass.to("g").magnitude

    result["time (h)"] = time

    result = pd.DataFrame.from_dict(result, orient='index').T

    os.makedirs(os.path.dirname(GROWTH_DATA_RAW), exist_ok=True)
    result.to_csv(GROWTH_DATA_CLEAN, index=False)


if __name__ == "__main__":
    main()
