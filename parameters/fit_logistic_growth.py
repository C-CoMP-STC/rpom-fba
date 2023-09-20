import json
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize

from parameters.drawdown import MASS_PER_CELL
from utils.math import get_interpolator, runge_kutta
from utils.units import u
from tqdm import tqdm

matplotlib.use("Agg")


UPTAKE_RATES_FILE = "data/drawdown_clean.csv"
LOGISTIC_RATES_OUTDIR = "parameters/logistic_fits/"
LOGISTIC_PLOTS_OUTDIR = "out/logistic_fits/"
RAW_UPTAKE_RATES_FILE = "data/drawdown_raw.xlsx"
RAW_UPTAKE_RATES_SHEET = "drawdown_clean"
GROWTH_DATA_FILE = "data/growth_curves_clean.csv"
GROWTH_RATES_FILE = "parameters/growth_rates/fitted_growth_rates.csv"
CUE_DATA_FILE = "data/CUE/cue_data.csv"
METHOD_NOT_RECOGNIZED = "{} is not a recognized method to fit uptake rates!"
BG_GRAY = "0.4"

import warnings
warnings.filterwarnings("error")

def logistic_growth_dynamic_system(S_0, X_0, k, alpha, max_t, dt=0.1):
    def dfdt(X):
        s, x = X
        try:
            return np.array([- alpha * k * s * x, k * s * x])
        except RuntimeWarning:
            pass

    time, trajectory, _ = runge_kutta(dfdt, np.array(
        [S_0, X_0]), 0, max_t, dt, pbar=False)

    return time, trajectory

def logistic_trajectory(t, S_0, X_0, k, alpha):
    r = k * S_0
    B = S_0 / alpha
    X_t = X_0 * B / (X_0 + (B - X_0) * np.exp(-r*t))
    S_t = alpha * (X_0 - X_t) + S_0

    return np.array([S_t, X_t])


def fit_logistic(t_S, S, t_X, X, dt=0.1, pbar=True):
    # def biomass_loss(params):
    #     k, alpha = params
    #     t, y_hat = logistic_growth_dynamic_system(
    #         S[0], X[0], k, alpha, t_S.max(), dt=dt)
    #     X_hat = get_interpolator(t, y_hat[:, 1])
    #     return np.sum((np.log(X_hat(t_X)) - np.log(X))**2)

    # def loss(params):
    #     k, alpha = params
    #     t, y_hat = logistic_growth_dynamic_system(
    #         S[0], X[0], k, alpha, t_S.max(), dt=dt)
    #     S_hat = get_interpolator(t, y_hat[:, 0])
    #     X_hat = get_interpolator(t, y_hat[:, 1])
    #     return np.sum((np.log(S_hat(t_S)) - np.log(S))**2) + np.sum((np.log(X_hat(t_X)) - np.log(X))**2)

    def loss(params):
        k, alpha = params
        S_hat = logistic_trajectory(t_S, S[0], X[0], k, alpha)[0, :]
        X_hat = logistic_trajectory(t_X, S[0], X[0], k, alpha)[1, :]
        return np.sum((np.log(S_hat) - np.log(S))**2) + np.sum((np.log(X_hat) - np.log(X))**2)

    with tqdm(desc="Fitting logistic: ") as pb:
        # k, alpha = minimize(biomass_loss,
        #                     [0.1, 1000],
        #                     callback=(lambda _: pb.update()) if pbar else None).x

        k, alpha = minimize(loss,
                            [0.1, 1000], # [k, alpha],
                            callback=(lambda _: pb.update()) if pbar else None).x
    
    return k, alpha


def plot_fits(
    mass_traces,
    substrate_traces,
    substrate_raw_data,
    rates
):
    # Ensure output directory exists
    os.makedirs(LOGISTIC_PLOTS_OUTDIR, exist_ok=True)

    for metabolite, (k, alpha) in rates.items():
        # Create plot
        fig, ax = plt.subplots()

        t_X, X = mass_traces[metabolite]
        t_S, S = substrate_traces[metabolite]

        # Plot points
        ax.scatter(*substrate_raw_data[metabolite], c=BG_GRAY)

        # Plot fitted curve
        t = np.linspace(0, t_S.max(), 100)

        # t, x = logistic_growth_dynamic_system(
        #     S[0],
        #     X[0],
        #     k,
        #     alpha,
        #     t_S.max())
        x = logistic_trajectory(t, S[0], X[0], k, alpha).T

        ax.plot(t, x[:, 0], "g--", alpha=0.5, label=f"[S]")
        ax.set_xlabel("Time (hr)")
        ax.set_ylabel("Substrate concentration (mM)")
        ax2 = ax.twinx()
        ax2.plot(t, x[:, 1], "r--")
        ax2.plot(t_X, X, color="r", label="X(t)")
        ax2.set_ylabel("Colony mass (g)", color="r")
        ax2.set_yscale("log")

        ax.set_title(
            f"{metabolite} \n($\\hat k = {k:.1f}\\ mM^{{-1}}h^{{-1}}$, $\\alpha={alpha:.2f} mM/gDCW$)")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("[Substrate] (mM)", color="g")

        # Save figure
        fig.set_size_inches(4, 2.25)
        fig.tight_layout()
        fig.savefig(os.path.join(LOGISTIC_PLOTS_OUTDIR,
                    f"uptake_[{metabolite}].png"))
        plt.close(fig)


def main():
    data = pd.read_csv(UPTAKE_RATES_FILE)
    cue_data = pd.read_csv(CUE_DATA_FILE)
    growth_data = pd.read_csv(GROWTH_DATA_FILE)

    # Filter cue_data down to just acetate singles
    acetate_data = (cue_data[(cue_data["Condition"] == "singles") &
                             (cue_data["Initial_mM_Acetate"] > 0)]
                    .groupby(["Time (h)", "Type"])["Value"]
                    .mean()
                    .reset_index())

    # Get mass traces (t, X)
    mass_traces = {
        metabolite: (growth_data["time (h)"][~growth_data[f"{metabolite}_predicted_mass"].isna()].values,
                     growth_data[f"{metabolite}_predicted_mass"].dropna().values)
        for metabolite in data["Compound"]
    }
    acetate_growth = acetate_data[acetate_data["Type"] == "counts"]
    mass_traces["acetate"] = (
        acetate_growth["Time (h)"].values, (acetate_growth["Value"].values * MASS_PER_CELL.to("g")).magnitude)

    # Get substrate traces (t, S)
    substrate_traces = {
        metabolite: (np.array([0, dt_hr]),
                     np.array([S_0, S_final]))
        for metabolite, dt_hr, S_0, S_final
        in zip(data["Compound"],
               data["dt_hr"],
               data["InitialMetabolite_mM"],
               data["FinalMetabolite_mM"])
    }
    substrate_traces["acetate"] = (acetate_data[acetate_data["Type"] == "drawdown (umol)"]["Time (h)"].values,
                                   (acetate_data[acetate_data["Type"] == "drawdown (umol)"]["Value"].values * u.umol).to("mmol").magnitude)

    rates = {
        metabolite: fit_logistic(substrate_traces[metabolite][0],
                                 substrate_traces[metabolite][1],
                                 t,
                                 X)
        for metabolite, (t, X)
        in mass_traces.items()
    }

    print("Plotting fitted rates...")
    raw_data = pd.read_excel(RAW_UPTAKE_RATES_FILE,
                             sheet_name=RAW_UPTAKE_RATES_SHEET)
    # Subset data to just one metabolite
    substrate_raw_data = {
        metabolite: (
            raw_data[raw_data["Compound"] == metabolite]["SampleType"].replace(
                {"initial": 0, "final": dt_hr}).values,
            mM_to_peak * raw_data[raw_data["Compound"]
                                  == metabolite]["IntegratedPeakArea"].values
        )
        for metabolite, dt_hr, mM_to_peak
        in zip(data["Compound"],
               data["dt_hr"],
               data["mM_to_peak_ratio"])
    }
    substrate_raw_data["acetate"] = (
        acetate_data[acetate_data["Type"] ==
                     "drawdown (umol)"]["Time (h)"].values,
        (acetate_data[acetate_data["Type"] == "drawdown (umol)"]
         ["Value"].values * u.umol).to("mmol").magnitude
    )

    plot_fits(mass_traces, substrate_traces, substrate_raw_data, rates)

    # Output rates to json
    print("Saving fitted rates to json...")
    os.makedirs(LOGISTIC_RATES_OUTDIR, exist_ok=True)

    # Writing to .json
    outfile = os.path.join(LOGISTIC_RATES_OUTDIR, "fitted_rates.json")
    with open(outfile, "w") as f:
        json.dump(rates, f, indent=4)


if __name__ == "__main__":
    main()
