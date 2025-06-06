import pickle
import numpy as np
import pandas as pd
from utils.units import u
from parameters.drawdown import MASS_PER_CELL, CUE_VOLUME


def main():
    CARBONS_PER_GLUCOSE = 6
    CARBONS_PER_ACETATE = 2
    OUTFILE = "data/clean/CUE/cue_data.csv"
    OUTFILE_DFBA = "data/clean/CUE/dFBA.pkl"

    files = {
        "doubles": {
            "counts": "data/raw/CUE/DSS3_doubles_counts.tsv",
            "cumulative": "data/raw/CUE/DSS3_doubles_cumulative.tsv",
            "drawdown": "data/raw/CUE/DSS3_doubles_drawdown.tsv",
            "flux": "data/raw/CUE/DSS3_doubles_flux.tsv",
            "OD": "data/raw/CUE/DSS3_doubles_OD600.tsv"
        },
        "singles": {
            "counts": "data/raw/CUE/DSS3_singles_counts.tsv",
            "cumulative": "data/raw/CUE/DSS3_singles_cumulative.tsv",
            "drawdown": "data/raw/CUE/DSS3_singles_drawdown.tsv",
            "flux": "data/raw/CUE/DSS3_singles_flux.tsv",
            "OD": "data/raw/CUE/DSS3_singles_OD600.tsv"
        }
    }

    data = {
        condition: {
            file: pd.read_csv(path, sep="\t")
            for file, path in file_dict.items()
        }
        for condition, file_dict in files.items()
    }

    # Clean counts (double)
    data["doubles"]["counts"] = data["doubles"]["counts"].melt(
        id_vars="Time (h)", var_name="Sample", value_name="Value")
    data["doubles"]["counts"]["Type"] = "counts"
    temp = data["doubles"]["counts"]["Sample"].str.extract(
        r"Rpom_T\d_(?P<Initial_C_Glucose>\dg)?(?P<Initial_C_Acetate>\da)?[nosub]?_?(?P<Sample>\d)?")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"4g": 4, "8g": 8, "4a": 4, "8a": 8, float("NaN"): 0})
    data["doubles"]["counts"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                              ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Clean counts (single)
    data["singles"]["counts"] = data["singles"]["counts"].melt(
        id_vars="Time (h)", var_name="Sample", value_name="Value")
    data["singles"]["counts"]["Type"] = "counts"
    temp = data["singles"]["counts"]["Sample"].str.extract(
        r"Rpom_(?P<Initial_C_Glucose>g)?(?P<Initial_C_Acetate>a)?(?P<Sample>\d)")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"g": 12, "a": 12, float("NaN"): 0})
    data["singles"]["counts"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                              ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Clean cumulative
    for condition in ["doubles", "singles"]:
        data[condition]["cumulative"] = data[condition]["cumulative"].melt(
            id_vars="Time (hours)", var_name="Sample", value_name="Value")
        data[condition]["cumulative"]["Type"] = "cumulative (umol)"
        data[condition]["cumulative"] = data[condition]["cumulative"].rename(
            columns={"Time (hours)": "Time (h)"})

    temp = data["doubles"]["cumulative"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>\dgluc)?_?(?P<Initial_C_Acetate>\dace)?[No_C_blank]?_?(?P<Sample>\d)?")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"4gluc": 4, "8gluc": 8, "4ace": 4, "8ace": 8, float("NaN"): 0})
    data["doubles"]["cumulative"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                                  ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    temp = data["singles"]["cumulative"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>Glc)?(?P<Initial_C_Acetate>Ac)?_(?P<Sample>[A-G])"
    )
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"Glc": 12, "Ac": 12, float("NaN"): 0})
    temp["Sample"] = temp["Sample"].replace(
        {"A": 1, "B": 2, "C": 3, "D": 1, "E": 2, "F": 3, "G": float("NaN")})
    data["singles"]["cumulative"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                                  ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Clean drawdown
    for condition in ["doubles", "singles"]:
        data[condition]["drawdown"].columns = data[condition]["drawdown"].columns.str.strip()
        data[condition]["drawdown"] = data[condition]["drawdown"].melt(
            id_vars="Time", var_name="Sample", value_name="Value")
        data[condition]["drawdown"]["Type"] = "drawdown (umol)"
        data[condition]["drawdown"] = data[condition]["drawdown"].rename(
            columns={"Time": "Time (h)"})

    temp = data["doubles"]["drawdown"]["Sample"].str.extract(
        r"R(?P<Initial_C_Glucose>\d)G(?P<Sample>\d)_(?P<Metabolite>acetate|glucose)")
    temp["Initial_C_Glucose"] = temp["Initial_C_Glucose"].astype(float)
    temp["Initial_C_Acetate"] = temp["Initial_C_Glucose"].replace({8: 4, 4: 8})
    data["doubles"]["drawdown"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample", "Metabolite"]
                                ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample", "Metabolite"]]

    temp = data["singles"]["drawdown"]["Sample"].str.extract(
        r"R[G|A](?P<Sample>\d)_(?P<Metabolite>acetate|glucose)")
    temp["Initial_C_Glucose"] = temp["Metabolite"].replace(
        {"acetate": 0., "glucose": 12.})
    temp["Initial_C_Acetate"] = temp["Metabolite"].replace(
        {"acetate": 12., "glucose": 0.})
    data["singles"]["drawdown"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample", "Metabolite"]
                                ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample", "Metabolite"]]

    # Clean flux
    for condition in ["doubles", "singles"]:
        data[condition]["flux"] = data[condition]["flux"].melt(
            id_vars="Time (hours)", var_name="Sample", value_name="Value")
        data[condition]["flux"]["Type"] = "flux (ppm/s)"
        data[condition]["flux"] = data[condition]["flux"].rename(
            columns={"Time (hours)": "Time (h)"})

    temp = data["doubles"]["flux"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>\d)gluc_(?P<Initial_C_Acetate>\d)ace_(?P<Sample>\d)")
    temp[["Initial_C_Glucose",
          "Initial_C_Acetate"]] = temp[["Initial_C_Glucose",
                                        "Initial_C_Acetate"]].astype("float").replace({float("NaN"): 0})
    data["doubles"]["flux"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                            ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    temp = data["singles"]["flux"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>Glc)?(?P<Initial_C_Acetate>Ac)?_(?P<Sample>[A-G])")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"Glc": 12.0, "Ac": 12.0, float("NaN"): 0.})
    temp["Sample"] = temp["Sample"].replace(
        {"A": 1, "B": 2, "C": 3, "D": 1, "E": 2, "F": 3, "G": float("NaN")})
    data["singles"]["flux"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                            ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Clean OD
    for condition in ["doubles", "singles"]:
        data[condition]["OD"] = data[condition]["OD"].melt(
            id_vars="Time", var_name="Sample", value_name="Value")
        data[condition]["OD"]["Type"] = "OD"
        data[condition]["OD"] = data[condition]["OD"].rename(
            columns={"Time": "Time (h)"})

    temp = data["doubles"]["OD"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>\d)gluc(?P<Initial_C_Acetate>\d)ace_(?P<Sample>\d)")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[[
        "Initial_C_Glucose", "Initial_C_Acetate"]].astype("float").replace({float("NaN"): 0})
    data["doubles"]["OD"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                          ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    temp = data["singles"]["OD"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>Glucose)?(?P<Initial_C_Acetate>Acetate)?_(?P<Sample>\d)")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"Glucose": 12.0, "Acetate": 12.0, float("NaN"): 0})
    data["singles"]["OD"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                          ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Merge and save all data in one file
    doubles = pd.concat(data["doubles"].values())
    doubles["Condition"] = "doubles"
    singles = pd.concat(data["singles"].values())
    singles["Condition"] = "singles"

    all_data = pd.concat([doubles, singles])
    all_data[["Initial_mM_Glucose", "Initial_mM_Acetate"]] = (
        all_data[["Initial_C_Glucose", "Initial_C_Acetate"]] / [CARBONS_PER_GLUCOSE, CARBONS_PER_ACETATE])
    all_data.to_csv(OUTFILE, index=False)

    # Create file with easily digestible condition data for running dFBA experiments
    dFBA_condition_data = {}
    initial_conditions = all_data[["Initial_mM_Glucose",
                                   "Initial_mM_Acetate"]].drop_duplicates().values

    for initial_glucose, initial_acetate in initial_conditions:
        condition_data = all_data[(all_data["Initial_mM_Glucose"] == initial_glucose) &
                                  (all_data["Initial_mM_Acetate"] == initial_acetate)]

        # Get biomass timeseries from data
        biomass_timeseries = condition_data[condition_data["Type"] == "counts"].groupby(
            "Time (h)")["Value"].mean().reset_index()
        biomass_timeseries["Value"] = (
            biomass_timeseries["Value"].values * (1/u.mL) * MASS_PER_CELL).to("g/L")
        b_t, b = biomass_timeseries.values.T

        initial_biomass = b[0]
        tmax = b_t.max()

        initial_glucose *= u.mM
        initial_acetate *= u.mM

        initial = np.array([
            initial_biomass,
            initial_glucose.magnitude,
            initial_acetate.magnitude
        ])

        # Get bounds from data
        substrate_timeseries = condition_data[condition_data["Type"] == "drawdown (umol)"].copy(
        )
        substrate_timeseries["Value"] = (
            (substrate_timeseries["Value"].values * u.umol) / CUE_VOLUME).to("mM")
        substrate_timeseries = substrate_timeseries.groupby(["Metabolite", "Time (h)"])[
            "Value"].mean().reset_index()

        g_t, g_s = substrate_timeseries[substrate_timeseries["Metabolite"] == "glucose"][[
            "Time (h)", "Value"]].values.T
        g_s /= CARBONS_PER_GLUCOSE
        a_t, a_s = substrate_timeseries[substrate_timeseries["Metabolite"] == "acetate"][[
            "Time (h)", "Value"]].values.T
        a_s /= CARBONS_PER_ACETATE
        if initial_glucose == 0:
            g_t, g_s = np.array([0, tmax]), np.array([0, 0])
        if initial_acetate == 0:
            a_t, a_s = np.array([0, tmax]), np.array([0, 0])

        # Raw biomass data
        mass_data = condition_data[condition_data["Type"] == "counts"].copy()
        mass_data["Mass (g/L)"] = ((mass_data["Value"].values /
                      u.mL) * MASS_PER_CELL).to("g/L")
        mass_data = mass_data.pivot(
            index="Time (h)",
            columns="Sample" if initial_glucose + initial_acetate != 0 else "Condition",
            values="Mass (g/L)").reset_index().values
        raw_b_t = mass_data[:, 0]
        raw_b = mass_data[:, 1:].T

        # Raw substrate data
        substrate_data = condition_data[condition_data["Type"]
                                        == "drawdown (umol)"].copy()
        substrate_data["Value"] = (
            (substrate_data["Value"].values * u.umol) / CUE_VOLUME).to("mM")

        raw_glucose = substrate_data[substrate_data["Metabolite"] == "glucose"].pivot(
            index="Time (h)", columns="Sample", values="Value").reset_index().values
        raw_g_t = raw_glucose[:,0]
        raw_g_s = raw_glucose[:,1:].T / CARBONS_PER_GLUCOSE

        raw_acetate = substrate_data[substrate_data["Metabolite"] == "acetate"].pivot(
            index="Time (h)", columns="Sample", values="Value").reset_index().values
        raw_a_t = raw_acetate[:,0]
        raw_a_s = raw_acetate[:,1:].T / CARBONS_PER_ACETATE

        dFBA_condition_data[(initial_glucose, initial_acetate)] = {
            "raw": {
                "raw_b_t": raw_b_t,
                "raw_b": raw_b,
                "raw_g_t": raw_g_t,
                "raw_g_s": raw_g_s,
                "raw_a_t": raw_a_t,
                "raw_a_s": raw_a_s,
            },
            "mean": {
                "g_t": g_t,
                "g_s": g_s,
                "a_t": a_t,
                "a_s": a_s,
                "b_t": b_t,
                "b_s": b
            }
        }
        
    # Save data for dFBA
    with open(OUTFILE_DFBA, "wb") as f:
        pickle.dump(dFBA_condition_data, f)
        

if __name__ == "__main__":
    main()
