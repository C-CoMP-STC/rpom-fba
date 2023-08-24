import matplotlib.pyplot as plt
import json
import os

import matplotlib
import numpy as np
import pandas as pd

from scipy.optimize import minimize

matplotlib.use("Agg")


UPTAKE_RATES_FILE = "data/drawdown_clean.csv"
UPTAKE_RATES_OUTDIR = "parameters/uptake_rates/"
UPTAKE_RATES_PLOTS_OUTDIR = "out/uptake_rates/"
RAW_UPTAKE_RATES_FILE = "data/drawdown_raw.xlsx"
RAW_UPTAKE_RATES_SHEET = "drawdown_clean"
GROWTH_DATA_FILE = "data/growth_curves_clean.csv"
GROWTH_RATES_FILE = "parameters/growth_rates/fitted_growth_rates.csv"
CUE_DATA_FILE = "data/CUE/cue_data.csv"
METHOD_NOT_RECOGNIZED = "{} is not a recognized method to fit uptake rates!"
BG_GRAY = "0.4"
MASS_PER_CELL = 0.95 # pg

# Michaelis-Menten
# Assuming K_M fixed for all metabolites, consider varying later
K_M = 1


def get_interpolator(t, y):
    def Y_t(times): return np.interp(np.atleast_1d(times), t, y)
    return Y_t


def michaelis_menten_dynamic_system(S_0, X_t, V_max, K_M, max_t, dt=0.1):
    time = np.arange(0, max_t, dt)
    trajectory = np.zeros((time.size, 2))

    s = S_0
    x = X_t(0)[0]
    for i, t in enumerate(time):
        trajectory[i, :] = [s, x]

        x = X_t(t)[0]
        dS = V_max * (s / (K_M + s)) * x

        s += dS * dt

    return np.array(time), np.array(trajectory)


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
        t, y_hat = michaelis_menten_dynamic_system(
            S[0], X_t, V_max, K_M, t_S.max(), dt=dt)
        S_hat = get_interpolator(t, y_hat[:, 0])
        return np.sum((S_hat(t_S) - S)**2)

    V_max = minimize(loss, [1]).x[0]
    return V_max


def fit_M_M(t_S, S, t_X, X, dt=0.1):
    """Fits a V_max and K_M, such that under Michaelis-Menten assumptions,
    the predicted drawdown curve of a substrate matches the given data.

    Args:
        t_S (array): Timepoints associated with S (must be same length as S).
        S (array): Substrate concentrations over time (must be same length as t_S).
        t_X (array): Timepoints associated with X (must be same length).
        X (array): Mass over time (must be same length as t_X).
        dt (float, optional): Timestep to use for the drawdown simulation. Defaults to 0.1.

    Returns:
        V_max: V_max that minimizes the squared sum loss.
        K_M: K_M that minimizes the squared sum loss.
    """

    X_t = get_interpolator(t_X, X)

    def loss(X):
        V_max, K_M = X
        t, y_hat = michaelis_menten_dynamic_system(
            S[0], X_t, V_max, K_M, t_S.max(), dt=dt)
        S_hat = get_interpolator(t, y_hat[:, 0])
        return np.sum((S_hat(t_S) - S)**2)

    V_max, K_M = minimize(loss, [1, 1]).x
    return V_max, K_M


def plot_fits(
    mass_traces,
    substrate_traces,
    substrate_raw_data,
    rates
):
    # Ensure output directory exists
    os.makedirs(UPTAKE_RATES_PLOTS_OUTDIR, exist_ok=True)

    for metabolite, rate in rates.items():
        # Create plot
        fig, ax = plt.subplots()

        t_X, X = mass_traces[metabolite]
        t_S, S = substrate_traces[metabolite]

        # Plot points
        ax.scatter(*substrate_raw_data[metabolite], c=BG_GRAY)

        # Plot fitted curve
        t = np.linspace(0, t_S.max(), 100)

        t, x = michaelis_menten_dynamic_system(
            S[0],
            get_interpolator(t_X, X),
            rate,
            K_M,
            t_S.max())

        ax.plot(t, x[:, 0], color="g", alpha=0.5, label=f"[S]")
        ax.set_xlabel("Time (hr)")
        ax.set_ylabel("Substrate concentration (mM)")
        ax2 = ax.twinx()
        ax2.plot(t, x[:, 1], color="r", label="X(t)")
        ax2.set_ylabel("Colony mass (g)", color="r")

        ax.set_title(
            f"{metabolite} \n($\\hat V_{{\\max}} = {abs(rate):.1f}\\ mmol / g \\cdot hr$)")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("[Substrate] (mM)", color="g")

        # Save figure
        fig.set_size_inches(4, 2.25)
        fig.tight_layout()
        fig.savefig(os.path.join(UPTAKE_RATES_PLOTS_OUTDIR,
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
        acetate_growth["Time (h)"].values, acetate_growth["Value"].values * MASS_PER_CELL * 1e-12)

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
                                   acetate_data[acetate_data["Type"] == "drawdown (umol)"]["Value"].values * 1e-3)

    rates = {
        metabolite: fit_M_M_for_given_K_M(substrate_traces[metabolite][0],
                                          substrate_traces[metabolite][1],
                                          t,
                                          X,
                                          K_M)
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
        acetate_data[acetate_data["Type"] == "drawdown (umol)"]["Time (h)"].values,
        acetate_data[acetate_data["Type"] == "drawdown (umol)"]["Value"].values * 1e-3
    )

    plot_fits(mass_traces, substrate_traces, substrate_raw_data, rates)

    # Output rates to json
    print("Saving fitted rates to json...")
    os.makedirs(UPTAKE_RATES_OUTDIR, exist_ok=True)

    # Writing to .json
    outfile = os.path.join(UPTAKE_RATES_OUTDIR, "fitted_uptake_rates.json")
    with open(outfile, "w") as f:
        json.dump(rates, f, indent=4)


if __name__ == "__main__":
    main()
