import pandas as pd


def main():
    CARBONS_PER_GLUCOSE = 6
    CARBONS_PER_ACETATE = 2

    files = {
        "doubles": {
            "counts": "data/CUE/DSS3_doubles_counts.tsv",
            "cumulative": "data/CUE/DSS3_doubles_cumulative.tsv",
            "drawdown": "data/CUE/DSS3_doubles_drawdown.tsv",
            "flux": "data/CUE/DSS3_doubles_flux.tsv",
            "OD": "data/CUE/DSS3_doubles_OD600.tsv"
        },
        "singles": {
            "counts": "data/CUE/DSS3_singles_counts.tsv",
            "cumulative": "data/CUE/DSS3_singles_cumulative.tsv",
            "drawdown": "data/CUE/DSS3_singles_drawdown.tsv",
            "flux": "data/CUE/DSS3_singles_flux.tsv",
            "OD": "data/CUE/DSS3_singles_OD600.tsv"
        }
    }

    data = {
        condition: {
            file: pd.read_csv(path, sep="\t")
            for file, path in file_dict.items()
        }
        for condition, file_dict in files.items()
    }

    # Clean counts
    data["doubles"]["counts"] = data["doubles"]["counts"].melt(
        id_vars="Time (h)", var_name="Sample", value_name="Value")
    data["doubles"]["counts"]["Type"] = "counts"
    temp = data["doubles"]["counts"]["Sample"].str.extract(
        r"Rpom_T\d_(?P<Initial_C_Glucose>\dg)?(?P<Initial_C_Acetate>\da)?[nosub]?_?(?P<Sample>\d)?")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"4g": 4, "8g": 8, "4a": 4, "8a": 8, float("NaN"): 0})
    data["doubles"]["counts"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                              ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Clean cumulative
    data["doubles"]["cumulative"] = data["doubles"]["cumulative"].melt(
        id_vars="Time (hours)", var_name="Sample", value_name="Value")
    data["doubles"]["cumulative"]["Type"] = "cumulative (umol)"
    data["doubles"]["cumulative"] = data["doubles"]["cumulative"].rename(
        columns={"Time (hours)": "Time (h)"})
    temp = data["doubles"]["cumulative"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>\dgluc)?_?(?P<Initial_C_Acetate>\dace)?[No_C_blank]?_?(?P<Sample>\d)?")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose", "Initial_C_Acetate"]].replace(
        {"4gluc": 4, "8gluc": 8, "4ace": 4, "8ace": 8, float("NaN"): 0})
    data["doubles"]["cumulative"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                                  ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Clean drawdown
    data["doubles"]["drawdown"] = data["doubles"]["drawdown"].melt(
        id_vars="Time", var_name="Sample", value_name="Value")
    data["doubles"]["drawdown"]["Type"] = "drawdown (umol)"
    data["doubles"]["drawdown"] = data["doubles"]["drawdown"].rename(
        columns={"Time": "Time (h)"})
    temp = data["doubles"]["drawdown"]["Sample"].str.extract(
        r"R(?P<Initial_C_Glucose>\d)G(?P<Sample>\d)_(?P<Metabolite>acetate|glucose)")
    temp["Initial_C_Glucose"] = temp["Initial_C_Glucose"].astype(float)
    temp["Initial_C_Acetate"] = temp["Initial_C_Glucose"].replace({8: 4, 4: 8})
    data["doubles"]["drawdown"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample", "Metabolite"]
                                ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample", "Metabolite"]]

    # Clean flux
    data["doubles"]["flux"] = data["doubles"]["flux"].melt(
        id_vars="Time (hours)", var_name="Sample", value_name="Value")
    data["doubles"]["flux"]["Type"] = "flux (ppm/s)"
    data["doubles"]["flux"] = data["doubles"]["flux"].rename(
        columns={"Time (hours)": "Time (h)"})
    temp = data["doubles"]["flux"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>\d)gluc_(?P<Initial_C_Acetate>\d)ace_(?P<Sample>\d)")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[["Initial_C_Glucose",
                                                             "Initial_C_Acetate"]].astype("float").replace({float("NaN"): 0})
    data["doubles"]["flux"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                            ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    # Clean OD
    data["doubles"]["OD"] = data["doubles"]["OD"].melt(
        id_vars="Time", var_name="Sample", value_name="Value")
    data["doubles"]["OD"]["Type"] = "OD"
    data["doubles"]["OD"] = data["doubles"]["OD"].rename(
        columns={"Time": "Time (h)"})
    temp = data["doubles"]["OD"]["Sample"].str.extract(
        r"(?P<Initial_C_Glucose>\d)gluc(?P<Initial_C_Acetate>\d)ace_(?P<Sample>\d)")
    temp[["Initial_C_Glucose", "Initial_C_Acetate"]] = temp[[
        "Initial_C_Glucose", "Initial_C_Acetate"]].astype("float").replace({float("NaN"): 0})
    data["doubles"]["OD"][["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]
                          ] = temp[["Initial_C_Glucose", "Initial_C_Acetate", "Sample"]]

    doubles = pd.concat(data["doubles"].values())
    doubles[["Initial_mM_Glucose", "Initial_mM_Acetate"]] = (
        doubles[["Initial_C_Glucose", "Initial_C_Acetate"]] / [CARBONS_PER_GLUCOSE, CARBONS_PER_ACETATE])
    doubles.to_csv("data/CUE/doubles.csv", index=False)


if __name__ == "__main__":
    main()
