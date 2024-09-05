import pickle
import numpy as np
import pandas as pd
from utils.units import u

from parameters.drawdown import (
    MASS_PER_CELL,
    MASS_PER_CELL_GLUCOSE,
    MASS_PER_CELL_ACETATE,
    CUE_VOLUME)

C_PER_GLUCOSE = 6
C_PER_ACETATE = 2
DATA_FILE = "data/clean/CUE2/cue2.xlsx"

cells = pd.read_excel(DATA_FILE, "Cell number", skiprows=1)
glucose = pd.read_excel(DATA_FILE, "Glucose", skiprows=1)
acetate = pd.read_excel(DATA_FILE, "Acetate", skiprows=1)
biomass_carbon = pd.read_excel(DATA_FILE, "C content", skiprows=1)
carbon = pd.read_excel(DATA_FILE, "CO2 cumulative", skiprows=1)
bge = pd.read_excel(DATA_FILE, "BGE", skiprows=1)


conditions = {
    "Glucose": (2.0 * u.mM, 0.0 * u.mM),
    "Acetate": (0.0 * u.mM, 6.0 * u.mM),
    "4gluc_8ace": ((4/C_PER_GLUCOSE) * u.mM, (8/C_PER_ACETATE) * u.mM),
    "8gluc_4ace": ((8/C_PER_GLUCOSE) * u.mM, (4/C_PER_ACETATE) * u.mM),
}

data = {}
for pref, key in conditions.items():
    g, a = key
    # Estimate the mass per cell by interpolating between glucose and acetate measurements
    # This is a very rough estimate, but it's the best we can do with the data we have
    theta = (a / (g + a)).magnitude
    mass_per_cell = np.interp(theta, [0, 1], [MASS_PER_CELL_GLUCOSE.magnitude, MASS_PER_CELL_ACETATE.magnitude]) * u.pg

    # Biomass
    cell_data = cells[[col for col in cells.columns if col.startswith(pref)]]  # Cells / mL
    b_t = cells["Time"].values
    cell_raw = cell_data[[
        col for col in cell_data.columns if not col.endswith("mean")]]
    cell_mean = cell_data[[
        col for col in cell_data.columns if col.endswith("mean")]]
    # Convert to g/L
    raw_b = (cell_raw.values * 1/u.mL * mass_per_cell).to("g/L").magnitude.T
    b_s = (cell_mean.values.T[0] * 1/u.mL * mass_per_cell).to("g/L").magnitude

    # Glucose
    gluc_data = glucose[[
        col for col in glucose.columns if col.startswith(pref)]]
    g_t = glucose["Time (single)" if pref in {
        "Glucose", "Acetate"} else "Time (double)"]
    g_t = g_t[~np.isnan(g_t)]
    glucose_raw = gluc_data.values[:g_t.size, :].T
    glucose_mean = glucose_raw.mean(axis=0)

    # Acetate
    ace_data = acetate[[
        col for col in acetate.columns if col.startswith(pref)]]
    a_t = acetate["Time (single)" if pref in {
        "Glucose", "Acetate"} else "Time (double)"]
    a_t = a_t[~np.isnan(a_t)]
    acetate_raw = ace_data.values[:a_t.size, :].T
    acetate_mean = acetate_raw.mean(axis=0)

    data[key] = {
        "raw": {
            "raw_b_t": b_t,
            "raw_b": raw_b,
            "raw_g_t": g_t,
            "raw_g_s": glucose_raw,
            "raw_a_t": a_t,
            "raw_a_s": acetate_raw
        },
        "mean": {
            "g_t": g_t,
            "g_s": glucose_mean,
            "a_t": a_t,
            "a_s": acetate_mean,
            "b_t": b_t,
            "b_s": b_s,
        }
    }

OUTFILE_DFBA = "data/clean/CUE2/dFBA.pkl"

# Save data for dFBA
with open(OUTFILE_DFBA, "wb") as f:
    pickle.dump(data, f)
