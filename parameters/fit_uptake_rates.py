import json
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


UPTAKE_RATES_FILE = "data/drawdown_clean.csv"
UPTAKE_RATES_OUTDIR = "parameters/uptake_rates/"
UPTAKE_RATES_PLOTS_OUTDIR = "out/uptake_rates/"
RAW_UPTAKE_RATES_FILE = "data/drawdown_raw.xlsx"
RAW_UPTAKE_RATES_SHEET = "drawdown_clean"
GROWTH_RATES_FILE = "parameters/growth_rates/fitted_growth_rates.csv"
METHOD_NOT_RECOGNIZED = "{} is not a recognized method to fit uptake rates!"
BG_GRAY = "0.4"

INITIAL_MASS = 1e-6  # Temporary guess, = 1ug


def fit_uptake_rates(
    data,
    growth_rates="GrowthRate",
    metabolites="Compound",
    initial_conc="InitialMetabolite_mM",
    final_conc="FinalMetabolite_mM",
    dts="dt_hr",
    method="linear"
):
    rates = {}

    # Prevent taking log of 0 using this method:
    def correct_zeros(x): return x if x != 0 else .00001

    match method:
        case "linear":
            rates = {
                metabolite: ((mu * (final - initial)) /
                             (INITIAL_MASS * (np.exp(mu * dt) - 1)))
                for metabolite, mu, final, initial, dt
                in zip(data[metabolites],
                       data[growth_rates],
                       data[final_conc],
                       data[initial_conc],
                       data[dts])
            }
        case "exponential":
            rates = {
                metabolite: (np.log(correct_zeros(cT)) -
                             np.log(correct_zeros(c0))) / dt
                for metabolite, c0, cT, dt
                in zip(
                    data[metabolites], +
                    data[initial_conc],
                    data[final_conc],
                    data[dts]
                )
            }
        case _:
            raise ValueError(METHOD_NOT_RECOGNIZED.format(method))

    return rates


def plot_fits(
    data,
    raw_data,
    rates,
    metabolites="Compound",
    dt="dt_hr",
    growth_rate="GrowthRate",
    method="linear"
):
    # Ensure output directory exists
    os.makedirs(UPTAKE_RATES_PLOTS_OUTDIR, exist_ok=True)

    for metabolite, rate in rates.items():
        # Subset data to just one metabolite
        metabolite_data = raw_data[raw_data[metabolites] == metabolite]

        total_time = data[data[metabolites] == metabolite][dt].values[0]
        mu = data[data[metabolites] == metabolite][growth_rate].values[0]
        mM_to_peak_ratio = data[data[metabolites] == metabolite]["mM_to_peak_ratio"].values[0]

        # Create plot
        fig, ax = plt.subplots()

        # Plot points
        t0 = metabolite_data[metabolite_data["SampleType"]
                             == "initial"]["IntegratedPeakArea"].values
        t_final = metabolite_data[metabolite_data["SampleType"]
                                  == "final"]["IntegratedPeakArea"].values
        t0 *= mM_to_peak_ratio
        t_final *= mM_to_peak_ratio
        ax.scatter(np.zeros_like(t0), t0, c=BG_GRAY)
        ax.scatter(np.ones_like(t_final) * total_time, t_final, c=BG_GRAY)

        mean_initial = t0.mean()

        # Plot fitted curve
        t = np.linspace(0, total_time, 100)
        match method:
            case "linear":
                ax.plot(t,
                        mean_initial + (rate * INITIAL_MASS * (np.exp(mu * t) - 1) / mu))
            case "exponential":
                ax.plot(t, mean_initial * np.exp(rate * t))
            case _:
                raise ValueError(METHOD_NOT_RECOGNIZED.format(method))

        ax.set_title(f"{metabolite} (fitted rate = {rate:.4f} $mmol / g \\cdot hr$)")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Substrate concentration (mM)")

        # Save figure
        fig.set_size_inches(8, 6)
        fig.tight_layout()
        fig.savefig(os.path.join(UPTAKE_RATES_PLOTS_OUTDIR,
                    f"uptake_[{metabolite}].png"))
        plt.close(fig)


def main():
    data = pd.read_csv(UPTAKE_RATES_FILE)
    growth_rates = (pd.read_csv(GROWTH_RATES_FILE)
                    .mean()
                    .reset_index(name="GrowthRate")
                    .rename(columns={"index": "Compound"}))
    data = data.merge(growth_rates, how="inner", on="Compound", validate="1:1")
    rates = fit_uptake_rates(data, method="linear")

    print("Plotting fitted rates...")
    raw_data = pd.read_excel(RAW_UPTAKE_RATES_FILE,
                             sheet_name=RAW_UPTAKE_RATES_SHEET)
    plot_fits(data, raw_data, rates)

    # Output rates to json
    print("Saving fitted rates to json...")
    os.makedirs(UPTAKE_RATES_OUTDIR, exist_ok=True)

    # Writing to .json
    outfile = os.path.join(UPTAKE_RATES_OUTDIR, "fitted_uptake_rates.json")
    with open(outfile, "w") as f:
        json.dump(rates, f, indent=4)


if __name__ == "__main__":
    main()
