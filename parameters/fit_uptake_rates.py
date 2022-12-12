import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

UPTAKE_RATES_FILE = "data/uptake_rates.xlsx"
UPTAKE_RATES_SHEET = "Eyeballed data"
UPTAKE_RATES_OUTDIR = "parameters/uptake_rates/"
UPTAKE_RATES_PLOTS_OUTDIR = "out/uptake_rates/"
METHOD_NOT_RECOGNIZED = "{} is not a recognized method to fit uptake rates!"
BG_GRAY = "0.4"


def load_data(filepath, sheet_name):
    data = pd.read_excel(filepath, sheet_name=sheet_name)
    return data


def fit_uptake_rates(
    data,
    metabolites="Metabolite",
    nmr_initial="T0_bar",
    nmr_final="WT_bar",
    initial_conc="Initial Concentration (mM metabolite)",
    dt="dt (hours)",
    method="exponential"
):
    rates = {}

    # Assuming NMR is linear with concentration,
    # convert final NMR value to concentration
    conc_to_nmr = data[initial_conc] / data[nmr_initial]
    conc_final = data[nmr_final] * conc_to_nmr

    # Prevent taking log of 0 using this method:
    def correct_zeros(x): return x if x != 0 else .00001

    match method:
        case "exponential":
            rates = {
                metabolite: (np.log(correct_zeros(cT)) -
                             np.log(correct_zeros(c0))) / dt
                for metabolite, c0, cT, dt
                in zip(
                    data[metabolites],
                    data[initial_conc],
                    conc_final,
                    data[dt]
                )
            }
        case _:
            raise ValueError(METHOD_NOT_RECOGNIZED.format(method))

    return rates


def plot_fits(
    data,
    rates,
    metabolites="Metabolite",
    nmr_initial_cols=["T0 (1)", "T0 (2)"],
    mean_nmr_initial="T0_bar",
    nmr_final_cols=["WT (1)", "WT (2)", "WT (3)"],
    mean_nmr_final="WT_bar",
    dt="dt (hours)",
    method="exponential"
):
    # Ensure output directory exists
    os.makedirs(UPTAKE_RATES_PLOTS_OUTDIR, exist_ok=True)

    for metabolite, rate in rates.items():
        # Subset data to just one metabolite
        metabolite_data = data[data[metabolites] == metabolite]

        # Extract useful variables
        mean_initial = metabolite_data[mean_nmr_initial].iloc[0]
        mean_final = metabolite_data[mean_nmr_final].iloc[0]
        total_time = metabolite_data[dt].iloc[0]

        # Create plot
        fig, ax = plt.subplots()

        # Plot points
        t0 = [metabolite_data[col] for col in nmr_initial_cols]
        t_final = [metabolite_data[col] for col in nmr_final_cols]
        ax.scatter(np.zeros_like(t0), t0, c=BG_GRAY)
        ax.scatter(np.ones_like(t_final) * total_time, t_final, c=BG_GRAY)

        # Plot mean of initial and final points as dotted horizontal lines
        ax.hlines([mean_initial, mean_final], 0, total_time, colors=[
                  BG_GRAY, BG_GRAY], linestyles="dashed")

        # Plot fitted curve
        t = np.linspace(0, total_time, 100)
        match method:
            case "exponential":
                ax.plot(t, mean_initial * np.exp(rate * t))
            case _:
                raise ValueError(METHOD_NOT_RECOGNIZED.format(method))

        ax.set_title(f"{metabolite} (fitted rate = {rate:.4f})")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Substrate concentration (area under NMR peak)")

        # Save figure
        fig.set_size_inches(8, 6)
        fig.tight_layout()
        fig.savefig(os.path.join(UPTAKE_RATES_PLOTS_OUTDIR,
                    f"uptake_[{metabolite}].png"))
        plt.close(fig)


def main():
    data = load_data(UPTAKE_RATES_FILE, UPTAKE_RATES_SHEET)
    rates = fit_uptake_rates(data)

    plot_fits(data, rates)

    # Output rates to csv
    os.makedirs(UPTAKE_RATES_OUTDIR, exist_ok=True)

    df = {"Metabolite": [], "Rate": []}
    for metabolite, rate in rates.items():
        df["Metabolite"].append(metabolite)
        df["Rate"].append(rate)
    df = pd.DataFrame(df)

    df.to_csv(
        os.path.join(UPTAKE_RATES_OUTDIR, "fitted_uptake_rates.csv"),
        index=False
    )


if __name__ == "__main__":
    main()
