import matplotlib.pyplot as plt
import json
import os

import matplotlib
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from tqdm import tqdm

from parameters.drawdown import K_M, MASS_PER_CELL, CUE_VOLUME, COLONY_VOLUME
from utils.math import get_interpolator, runge_kutta
from utils.units import u

matplotlib.use("Agg")


def michaelis_menten_dynamic_system(S_0, X_t, V_max, K_M, max_t, dt=0.1):
    def df_dt(y):
        s, t = y
        dS = V_max * (s / (K_M + s)) * X_t(t)[0]

        # Slightly hacky, using time as a state variable with dT/dt=1
        return np.array([dS, 1])

    t, trajectory, _ = runge_kutta(
        df_dt, np.array([S_0, 0]), 0, max_t, dt, pbar=False)
    S_t = trajectory[:, 0][..., np.newaxis]
    return t, np.concatenate([S_t, X_t(t)[..., np.newaxis]], axis=1)


def fit_M_M_for_given_K_M(t_S, S, t_X, X, K_M, dt=0.1):
    """Fits a V_max for the given K_M, such that under Michaelis-Menten assumptions,
    the predicted drawdown curve of a substrate matches the given data.

    Args:
        t_S (array): Timepoints associated with S (must be same length as S).
        S (array): Substrate concentrations over time (must be same length as t_S).
        t_X (array): Timepoints associated with X (must be same length).
        X (array): Mass over time (must be same length as t_X).
        K_M (float): Fixed K_M over which to fit V_max.
        dt (float, optional): Timestep to use for the drawdown simulation. Defaults to 0.1.

    Returns:
        V_max: V_max that minimizes the squared sum loss.
    """

    X_t = get_interpolator(t_X, X)

    def loss(V_max):
        V_max = V_max[0]
        t, y_hat = michaelis_menten_dynamic_system(
            S[0], X_t, V_max, K_M, t_S.max(), dt=dt)
        S_hat = get_interpolator(t, y_hat[:, 0])
        return np.sum((S_hat(t_S) - S)**2)

    V_max = minimize(loss, [1]).x[0]
    return V_max


def plot_fits(
    metabolite_traces,
    substrate_raw_data,
    rates,
    out_dir
):
    for metabolite, rate in rates.items():
        traces = metabolite_traces[metabolite]
        t_X, X = traces["t_X"], traces["X"]
        t_S, S = traces["t_S"], traces["S"]
        data_t, data_S = substrate_raw_data[metabolite]

        fig, _ = plot_fit(metabolite, rate, t_X, X, t_S, S,
                          data_t.magnitude, data_S.magnitude)
        fig.set_size_inches(4, 2.25)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"uptake_[{metabolite}].png"))
        plt.close(fig)


def plot_fit(metabolite, rate, t_X, X, t_S, S, data_t, data_S):
    BG_GRAY = "0.4"

    fig, ax = plt.subplots()

    # Plot fitted curve
    t, x = michaelis_menten_dynamic_system(
        S[0].magnitude,
        get_interpolator(t_X.magnitude, X.magnitude),
        rate,
        K_M.magnitude,
        t_S.max().magnitude)

    # Plot datapoints, simulation
    ax.scatter(data_t, data_S, c=BG_GRAY)
    ax.plot(t, x[:, 0], color="g", alpha=0.5, label=f"[S]")

    # Plot biomass data
    ax2 = ax.twinx()
    ax2.plot(t, x[:, 1], color="r", label="X(t)")

    # Labels
    ax.set_title(
        f"{metabolite} \n($\\hat V_{{\\max}} = {abs(rate):.2f}\\ mmol/g\\cdot hr$)")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("[Substrate] (mM)", color="g")
    ax2.set_ylabel("Colony mass (g/L)", color="r")

    return fig, ax


def main():
    UPTAKE_RATES_FILE = "data/drawdown_clean.csv"
    UPTAKE_RATES_OUTDIR = "parameters/uptake_rates/"
    UPTAKE_RATES_PLOTS_OUTDIR = "out/uptake_rates/"
    RAW_UPTAKE_RATES_FILE = "data/drawdown_raw.xlsx"
    RAW_UPTAKE_RATES_SHEET = "drawdown_clean"
    GROWTH_DATA_FILE = "data/growth_curves_clean.csv"
    CUE_DATA_FILE = "data/CUE/cue_data.csv"

    drawdown_data = pd.read_csv(UPTAKE_RATES_FILE)
    cue_data = pd.read_csv(CUE_DATA_FILE)
    growth_data = pd.read_csv(GROWTH_DATA_FILE)
    raw_data = pd.read_excel(RAW_UPTAKE_RATES_FILE,
                             sheet_name=RAW_UPTAKE_RATES_SHEET)

    # Filter cue_data down to just acetate singles
    acetate_data = (cue_data[(cue_data["Condition"] == "singles") &
                             (cue_data["Initial_mM_Acetate"] > 0)]
                    .groupby(["Time (h)", "Type"])["Value"]
                    .mean()
                    .reset_index())
    acetate_growth = acetate_data[acetate_data["Type"] == "counts"]
    acetate_growth["Value"] = ((acetate_growth["Value"].values / u.mL ) * CUE_VOLUME).to("dimensionless").magnitude

    # Get data for each metabolite, namely biomass traces (t_X, X) and substrate traces (t_S, S)
    metabolite_traces = {
        metabolite: {
            "t_X": growth_data["time (h)"][~growth_data[f"{metabolite}_predicted_mass"].isna()].values * u.hr,
            "X": growth_data[f"{metabolite}_predicted_mass"].dropna().values * u.g / COLONY_VOLUME.to("L"),
            "t_S": np.array([0, dt_hr]) * u.hr,
            "S": np.array([S_0, S_final]) * u.mM
        }
        for metabolite, dt_hr, S_0, S_final in zip(drawdown_data["Compound"],
                                                   drawdown_data["dt_hr"],
                                                   drawdown_data["InitialMetabolite_mM"],
                                                   drawdown_data["FinalMetabolite_mM"])
    }
    metabolite_traces["acetate"] = {
        "t_X": acetate_growth["Time (h)"].values * u.hr,
        "X": (acetate_growth["Value"].values * MASS_PER_CELL).to("g") / CUE_VOLUME.to("L"),
        "t_S": acetate_data[acetate_data["Type"] == "drawdown (umol)"]["Time (h)"].values * u.hr,
        "S": (acetate_data[acetate_data["Type"] == "drawdown (umol)"]["Value"].values * u.umol / CUE_VOLUME).to("mM")
    }

    # Run fitting
    rates = {
        metabolite: fit_M_M_for_given_K_M(trace["t_S"].magnitude,
                                          trace["S"].magnitude,
                                          trace["t_X"].magnitude,
                                          trace["X"].magnitude,
                                          K_M.magnitude)
        for metabolite, trace
        in tqdm(metabolite_traces.items(), desc="Fitting")
    }

    # Get raw points for plot
    substrate_raw_data = {
        metabolite: (
            raw_data[raw_data["Compound"] == metabolite]["SampleType"].replace(
                {"initial": 0, "final": dt_hr}).values * u.hr,
            mM_to_peak * raw_data[raw_data["Compound"]
                                  == metabolite]["IntegratedPeakArea"].values * u.mM
        )
        for metabolite, dt_hr, mM_to_peak
        in tqdm(zip(drawdown_data["Compound"],
                    drawdown_data["dt_hr"],
                    drawdown_data["mM_to_peak_ratio"]),
                desc="Plotting")
    }
    substrate_raw_data["acetate"] = (
        acetate_data[acetate_data["Type"] ==
                     "drawdown (umol)"]["Time (h)"].values * u.hr,
        (acetate_data[acetate_data["Type"] ==
                      "drawdown (umol)"]["Value"].values * u.umol / CUE_VOLUME).to("mM")
    )

    # Make and save plots
    os.makedirs(UPTAKE_RATES_PLOTS_OUTDIR, exist_ok=True)
    plot_fits(metabolite_traces, substrate_raw_data,
              rates, UPTAKE_RATES_PLOTS_OUTDIR)

    # Output rates to json
    os.makedirs(UPTAKE_RATES_OUTDIR, exist_ok=True)
    with open(os.path.join(UPTAKE_RATES_OUTDIR, "fitted_uptake_rates.json"), "w") as f:
        json.dump(rates, f, indent=4)


if __name__ == "__main__":
    main()
