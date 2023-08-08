"""
TODO: Maybe scrap this, approach I was thinking of taking shouldn't work
"""


import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from parameters.od_calibration import get_od_to_cell_count_calibration


def main():
    GROWTH_CURVES_FILE = "data/growth_curves_raw.csv"
    OUTDIR = "out/maintenance"
    TIME = "time (h)"

    # Load data
    data = pd.read_csv(GROWTH_CURVES_FILE, skiprows=[1])

    # Ensure output directory exists
    os.makedirs(OUTDIR, exist_ok=True)

    # Convert from OD to cell counts
    od_to_count = get_od_to_cell_count_calibration()
    data.loc[:, data.columns != TIME] = (data.loc[:, data.columns != TIME].applymap(
        lambda s: od_to_count.predict(s).flatten()[0], na_action="ignore"))

    compounds = {s[:(s.index(".") if "." in s else None)]
                 for s in data.columns if s not in {TIME, "mean blank"}}

    for compound in compounds:
        # Get data for specified compound
        compound_data = data[[
            TIME, *[c for c in data.columns if c.startswith(compound)]]]
        compound_data = compound_data.dropna(subset=[c for c in compound_data.columns if c != time])

        # Get estimates of growth rate mu from each succesive pair of timepoints
        dts = compound_data[TIME].diff()
        mu_hats = np.log(compound_data.loc[:, compound_data.columns != TIME]).diff(
            axis=0).apply(lambda c: c / dts)


if __name__ == "__main__":
    main()
