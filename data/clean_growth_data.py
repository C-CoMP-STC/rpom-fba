import pickle
import numpy as np
import pandas as pd

from parameters.drawdown import MASS_PER_CELL

GROWTH_DATA_FILE = "data/growth_curves_raw.csv"
OUT_FILE = "data/growth_curves_clean.csv"

def main():
    # Need to:
    # - convert metabolite names to recognizable IDs
    # - convert OD to predicted cell count and predicted mass
    data = pd.read_csv(GROWTH_DATA_FILE)

    with open("parameters/conversions/od_to_cell_count.pickle", "rb") as f:
        od_reg = pickle.load(f)

    carbon_sources = {
        'glycerol',
        'carnitine',
        'malate',
        'fumarate',
        'succinate',
        'citrate',
        '3-OH_butyrate',
        'ectoine',
        'DMSP',
        'thymidine',
        'cysteate',
        'xylose',
        'glucose',
        'isethionate',
        'cadaverine',
        'putresceine',
        'spermidine',
        'DHPS',
        'taurine',
        'choline',
        'GlcNac'
    }

    blank = data["mean blank"].values[1:].astype(float)
    time = data["time (h)"].values[1:].astype(float)

    result = {}
    for carbon_source in carbon_sources:
        # Subset to data pertaining to this carbon source
        source_data = data.loc[:, [c.startswith(carbon_source) for c in data.columns]].values[1:].astype(float)

        # Take the mean across rows
        mean_growth = source_data.mean(axis=1)
        result[f"{carbon_source}_mean_OD"] = mean_growth

        # Translate to cell count
        cell_count = od_reg.predict(mean_growth[~np.isnan(mean_growth)]).T[0]
        result[f"{carbon_source}_predicted_count"] = cell_count

        # Translate to mass
        total_mass = MASS_PER_CELL * cell_count
        result[f"{carbon_source}_predicted_mass"] = total_mass / 1e12  # convert to g

    result["time (h)"] = time

    result = pd.DataFrame.from_dict(result, orient='index').T
    result.to_csv(OUT_FILE, index=False)


if __name__ == "__main__":
    main()
