import matplotlib.pyplot as plt
import os

import matplotlib
import numpy as np
import pandas as pd

from parameters.od_calibration import get_od_to_cell_density_calibration, CellDensityRegressor
from parameters.drawdown import COLONY_VOLUME
from utils.units import u

matplotlib.use("Agg")

GROWTH_CURVES_FILE = "data/growth_curves_raw.csv"
GROWTH_RATE_OUTDIR = "parameters/growth_rates"
GROWTH_RATE_PLOT_OUTDIR = "out/growth_rates"


def fit_growth_curve(data, compound, time="time (h)"):
    # Discard any rows from after end of experiment
    data = data.dropna(subset=[c for c in data.columns if c != time])

    # Get estimates of growth rate mu from each succesive pair of timepoints
    dts = data[time].diff()
    mu_hats = np.log(data.loc[:, data.columns != time]).diff(
        axis=0).apply(lambda c: c / dts)

    # Take median as the final estimate for each replicate
    median_mus = mu_hats.median()

    # Plot distibution of mu_hats over time for each replicate
    fig, axs = plt.subplots(len(data.columns) - 1, 1)
    for ax, replicate, estimate in zip(axs, mu_hats.values.T, median_mus):
        ax2 = ax.twiny()

        # Plot histogram of mu_hats
        ax2.hist(replicate, alpha=0.25, orientation="horizontal")
        ax2.set_xticks([])

        # Plot mu vs time
        ax.plot(data[time], replicate)
        ax.set_ylabel(r"$\hat\mu$")

        # Plot estimated mu
        ax.hlines(
            y=estimate,
            xmin=data[time].min(),
            xmax=data[time].max(),
            colors="k",
            linestyles="dashed",
            label=f"Median $\\mu={estimate:.3f}$"
        )

        ax.legend()

    axs[0].set_title(compound)
    axs[-1].set_xlabel("Time (h)")

    return median_mus, fig, axs


def main():
    # Load data
    data = pd.read_csv(GROWTH_CURVES_FILE, skiprows=[1])

    # Convert from OD to cell counts
    od_to_cell_density = get_od_to_cell_density_calibration()
    data.loc[:, data.columns != "time (h)"] = (
        data.loc[:, data.columns != "time (h)"].applymap(
            lambda s: ((od_to_cell_density.predict(s).flatten()[
                       0] / u.mL) * COLONY_VOLUME).to("dimensionless").magnitude,
            na_action="ignore")
    )

    compounds = {s[:(s.index(".") if "." in s else None)]
                 for s in data.columns if s not in {"time (h)", "mean blank"}}

    os.makedirs(GROWTH_RATE_PLOT_OUTDIR, exist_ok=True)

    growth_rates = {}
    for compound in compounds:
        estimates, fig, _ = fit_growth_curve(
            data[["time (h)", *[c for c in data.columns if c.startswith(compound)]]], compound)

        # Store estimated growth rates
        growth_rates[compound] = estimates.values

        # Save plot
        fig.set_size_inches(8, 6)
        fig.tight_layout()
        fig.savefig(os.path.join(GROWTH_RATE_PLOT_OUTDIR, f"{compound}.png"))
        plt.close(fig)

    # Save estimates
    os.makedirs(GROWTH_RATE_OUTDIR, exist_ok=True)
    pd.DataFrame(growth_rates).to_csv(
        os.path.join(GROWTH_RATE_OUTDIR, "fitted_growth_rates.csv"),
        index=False)


if __name__ == "__main__":
    main()
