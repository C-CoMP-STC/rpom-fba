import pandas as pd

from data.cleaning.files import DRAWDOWN_DATA_FILE, UPTAKE_RATES_SHEET, METADATA_SHEET

def clean_drawdown_data(filepath):
    data = pd.read_excel(filepath, sheet_name=UPTAKE_RATES_SHEET)
    meta = pd.read_excel(filepath, sheet_name=METADATA_SHEET)

    # Calculate conversion factors for each metabolite between
    # concentration (mM) and integrated peak area
    # 
    # NOTE: Assumes linear relationship between NMR peak area and concentration
    result = (
        data[data["SampleType"] == "initial"]
        .groupby("Compound")["IntegratedPeakArea"]
        .mean()
        .reset_index()
        .merge(
            meta[["Metabolite", "InitialMetabolite_mM", "dt_hr"]],
            how="inner",
            left_on="Compound",
            right_on="Metabolite",
            validate="1:1")
        .drop(columns="Metabolite")
    )
    result["mM_to_peak_ratio"] = (
        result["InitialMetabolite_mM"] / result["IntegratedPeakArea"])

    # Calculate final concentrations in mM
    result = result.merge(
        data[data["SampleType"] == "final"]
        .groupby("Compound")["IntegratedPeakArea"]
        .mean()
        .reset_index(),
        how="inner",
        on="Compound",
        suffixes=["_initial", "_final"],
        validate="1:1")
    result["FinalMetabolite_mM"] = result["IntegratedPeakArea_final"] * result["mM_to_peak_ratio"]

    return result


def main():
    OUTFILE = "data/drawdown_clean.csv"
    clean_data = clean_drawdown_data(DRAWDOWN_DATA_FILE)
    clean_data.to_csv(OUTFILE, index=False)


if __name__ == "__main__":
    main()
