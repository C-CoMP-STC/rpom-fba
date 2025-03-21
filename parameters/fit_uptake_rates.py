import json
import os
import pickle

import matplotlib
import numpy as np
import pandas as pd
from scipy.optimize import minimize, fmin
from scipy.integrate import odeint
from tqdm import tqdm

from data.files import (DRAWDOWN_DATA_CLEAN, DRAWDOWN_DATA_RAW,
                        GROWTH_DATA_CLEAN, UPTAKE_RATES_SHEET)
from parameters.drawdown import COLONY_VOLUME, CUE_VOLUME, K_M, MASS_PER_CELL
from utils.colors import *
from utils.math import get_interpolator, runge_kutta
from utils.units import u

matplotlib.use("Agg")

# To fit the uptake parameters, we first need to extract the exponential phase of growth,
# where the uptake rate is likely saturated. 
# Below, I use a method where I look for regions where the mean instantaneous growth rate is greater than a threshold,
# and the coefficient of variation (CV) of instantaneous growth rates across replicates is less than a threshold.
NO_GROWTH_THRESHOLD = 0.01
CV_THRESHOLD = 0.35

DETECTION_LIMIT = 1e-4  # Constant to add to prevent log of zero


def plot_data(ax, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates, lw=2, N_color = C_GLUCOSE):
    ax.plot(t_X_obs, X_obs_replicates.T, marker="o", lw=0, color=C_BIOMASS, alpha=0.25)
    ax.plot(t_X_obs, X_obs_replicates.mean(axis=0), color=C_BIOMASS)
    ax2 = ax.twinx()
    ax2.plot(t_N_obs, N_obs_replicates.T, marker="o", lw=lw, color=N_color)

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Biomass (g/L)")
    ax2.set_ylabel("Substrate (mM)")

    return ax, ax2

def get_exponential_region(t_X_obs, X_obs_replicates, CV_THRESHOLD=CV_THRESHOLD, NO_GROWTH_THRESHOLD=NO_GROWTH_THRESHOLD):
    # find mu = dlog(X)/dt
    X_obs_replicates = np.log(X_obs_replicates + DETECTION_LIMIT)
    mean_dX_dt = np.mean([np.gradient(X_obs, t_X_obs) for X_obs in X_obs_replicates], axis=0)

    # Get the standard deviation, CV of the gradient
    dX_dt_std = np.std([np.gradient(X_obs, t_X_obs) for X_obs in X_obs_replicates], axis=0)
    dX_dt_CV = dX_dt_std / np.abs(mean_dX_dt)

    # Extract the exponential region according to CV and mu thresholds
    exponential_region = (dX_dt_CV < CV_THRESHOLD) & (mean_dX_dt > NO_GROWTH_THRESHOLD)

    # Filter to longest continuous region
    region_starts = []
    region_lengths = []
    in_region = False
    for i, exp in enumerate(exponential_region):
        if exp and not in_region:
            in_region = True
            region_starts.append(i)
        elif not exp and in_region:
            in_region = False
            region_lengths.append(i - region_starts[-1])

    longest = np.argmax(region_lengths)
    exponential_region = np.zeros_like(exponential_region, dtype=bool)
    exponential_region[region_starts[longest] : region_starts[longest] + region_lengths[longest]] = True
    
    return exponential_region


def plot_exponential_region_finder(ax, t_X_obs, X_obs_replicates, exponential_region, CV_THRESHOLD=CV_THRESHOLD, NO_GROWTH_THRESHOLD=NO_GROWTH_THRESHOLD):
    # Log-transform
    X_obs_replicates = np.log(X_obs_replicates + DETECTION_LIMIT)

    # Plot mu datapoints
    for X_obs in X_obs_replicates:
        dX_dt = np.gradient(X_obs, t_X_obs)
        ax.hlines(0, 0, t_X_obs[-1], color="0.5", linestyle="-", zorder=-1)
        ax.plot(t_X_obs, dX_dt, marker="o", lw=0, color=C_BIOMASS, alpha=0.25)
    
    # Plot mean
    mean_dX_dt = np.mean([np.gradient(X_obs, t_X_obs) for X_obs in X_obs_replicates], axis=0)
    ax.plot(t_X_obs, mean_dX_dt, color=C_BIOMASS)

    # Get the standard deviation of the gradient
    ax2 = ax.twinx()
    dX_dt_std = np.std([np.gradient(X_obs, t_X_obs) for X_obs in X_obs_replicates], axis=0)
    dX_dt_CV = dX_dt_std / np.abs(mean_dX_dt)
    ax2.plot(t_X_obs, dX_dt_CV, color="r")
    ax2.set_ylim(0, 1)
    ax2.set_ylabel(r"$\sigma/\mu$", color="r")

    # Plot the exponential region
    t_exponential = t_X_obs[exponential_region]
    growth_rate = np.mean(mean_dX_dt[exponential_region])
    ax.fill_between(t_X_obs, -1, 1, where=exponential_region, color="0.9", zorder=-2)
    ax.hlines(growth_rate, t_exponential[0], t_exponential[-1], color="k", linestyle="--", zorder=-1)
    ax.text(t_exponential[0] + 2, growth_rate + 0.1, f"μ = {growth_rate:.2f}", horizontalalignment="left")

    # Plot thresholds
    ax.hlines(NO_GROWTH_THRESHOLD, 0, t_X_obs[-1], color=C_BIOMASS, linestyle="--", zorder=-1)
    ax2.hlines(CV_THRESHOLD, 0, t_X_obs[-1], color="r", linestyle="--", zorder=-1)

    ax.set_xlim(0, t_X_obs[-1])
    ax.set_ylim(-0.1, 0.8)
    ax.set_ylabel(r"$\frac{{d}}{{dt}}\log (X)$", color=C_BIOMASS)
    ax.set_xlabel("Time (h)")


def fit_exp_region_params(t_N_obs, N_obs_replicates, t_X_obs, X_obs_replicates, exponential_region):
    # Get mean growth rate in exponential region (mu)
    mean_dX_dt = np.mean([np.gradient(np.log(X_obs), t_X_obs) for X_obs in X_obs_replicates], axis=0)
    mu = np.mean(mean_dX_dt[exponential_region])
    
    # Filter to exponential region
    t_exponential = t_X_obs[exponential_region]
    X_obs_replicates = X_obs_replicates[:, exponential_region]
    X0 = X_obs_replicates.mean(axis=0)[0]
    t_shifted = t_exponential - t_exponential[0]

    # Make nutrient interpolator
    N_interp = []
    for N_obs in N_obs_replicates:
        interp = get_interpolator(t_N_obs, N_obs)
        N_interp.append(interp(t_X_obs))
    N_interp = np.array(N_interp)
    N_interp = N_interp[:, exponential_region]

    # def get_model_error(t, N, X0, mu):
    #     def model_error(N0, k):
    #         return np.sum((N - (N0 - k * X0 * (np.exp(mu * t) - 1) / mu))**2)
    #     return model_error

    def model_error(N0, k):
        N_true = N_interp.mean(axis=0)
        return np.sum((N_true - (N0 - k * X0 * (np.exp(mu * t_shifted) - 1) / mu))**2)

    # Fit N0, Vmax
    # err = get_model_error(t_shifted, N_interp.mean(axis=0), X0, mu)
    sol = minimize(lambda x: model_error(*x), [N_interp.mean(axis=0)[0], 1])
    N0, Vmax = sol.x

    return mu, N0, Vmax


def plot_exp_region_params(ax, mu, N0, Vmax, t_N_obs, N_obs_replicates, t_X_obs, X_obs_replicates, exponential_region, nutrient_color = C_GLUCOSE):
    # Make nutrient interpolator
    N_interp = []
    for N_obs in N_obs_replicates:
        interp = get_interpolator(t_N_obs, N_obs)
        N_interp.append(interp(t_X_obs))
    N_interp = np.array(N_interp)
    
    # Filter to exponential region
    t_exponential = t_X_obs[exponential_region]
    t_shifted = t_exponential - t_exponential[0]
    t_N_exponential = t_N_obs[((t_exponential[0] < t_N_obs) & (t_N_obs < t_exponential[-1])).values]
    X_obs_replicates = X_obs_replicates[:, exponential_region]
    N_obs_replicates = N_obs_replicates[:, ((t_exponential[0] < t_N_obs) & (t_N_obs < t_exponential[-1])).values]
    N_interp = N_interp[:, exponential_region]

    # Plot data    
    _, ax2 = plot_data(ax, t_exponential, X_obs_replicates, t_N_exponential, N_obs_replicates, N_color = nutrient_color, lw=0)
    ax2.plot(t_exponential, N_interp.T, color=nutrient_color, lw=1, alpha=0.5, marker="x")

    # Plot the model
    X0 = X_obs_replicates.mean(axis=0)[0]
    model_fit = N0 - Vmax * X0 * (np.exp(mu * t_shifted) - 1) / mu
    ax.plot(t_exponential, X0 * np.exp(mu * t_shifted), color="k", linestyle="--")
    ax2.plot(t_exponential, model_fit, color="k", linestyle="--")
    ax2.text(t_exponential[0] + 2, N0, f"V_max = {Vmax:.2f}", horizontalalignment="left")

    return ax, ax2


def plot_error_vs_K_M(ax, V_max, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates):
    # Take mean of X
    X_obs_mean = X_obs_replicates.mean(axis=0)

    # Define derivative
    def michaelis_menten(N, t, K_M, X=X_obs_mean, t_X_obs=t_X_obs):
        # equations
        X_t = np.interp(t, t_X_obs, X)
        return -V_max * (N / (K_M + N)) * X_t
    
    # Define error function
    def error(log_K_M):
        K_M = np.exp(log_K_M)
        N_model = odeint(michaelis_menten, N_obs_replicates.mean(axis=0)[0], t_N_obs, args=(K_M,)).flat
        # TODO: don't average before taking error?
        return np.sum((N_obs_replicates.mean(axis=0) - N_model)**2)

    # Calculate errors
    K_Ms = np.linspace(0.001, 1, 100)
    log_K_Ms = np.log(K_Ms)
    errors = np.array([error(log_K_M) for log_K_M in log_K_Ms])

    # Plot
    ax.plot(K_Ms, errors)
    ax.set_xlabel("K_M")
    ax.set_ylabel("Error")

    return ax


def fit_K_M(V_max, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates):
    X_obs_mean = X_obs_replicates.mean(axis=0)

    # Define derivative
    def michaelis_menten(N, t, K_M, X=X_obs_mean, t_X_obs=t_X_obs):
        # equations
        X_t = np.interp(t, t_X_obs, X)
        return -V_max * (N / (K_M + N)) * X_t
    
    # Define error function
    def error(log_K_M):
        K_M = np.exp(log_K_M)
        N_model = odeint(michaelis_menten, N_obs_replicates.mean(axis=0)[0], t_N_obs, args=(K_M,)).flat
        # TODO: don't average before taking error?
        return np.sum((N_obs_replicates.mean(axis=0) - N_model)**2)

    # Fit the model
    log_K_M = fmin(error, np.log(0.1))[0]
    K_M = np.exp(log_K_M)

    return K_M


def plot_K_M_fit(ax, V_max, K_M, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates, nutrient_color = C_GLUCOSE):
    X_obs_mean = X_obs_replicates.mean(axis=0)
    
    # Define derivative
    def michaelis_menten(N, t, K_M, X=X_obs_mean, t_X_obs=t_X_obs):
        # equations
        X_t = np.interp(t, t_X_obs, X)
        return -V_max * (N / (K_M + N)) * X_t
    
    # Plot raw data
    _, ax2 = plot_data(ax, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates, N_color = nutrient_color, lw=0)

    # Plot the model
    N_model = odeint(michaelis_menten, N_obs_replicates.mean(axis=0)[0], t_X_obs, args=(K_M,))
    ax2.plot(t_X_obs, N_model, color="k", linestyle="--")
    ax2.text(2, 0.9 * N_model[0], f"$V_{{max}} = {V_max:.2f}$\n$K_M = {K_M:.2g}$", horizontalalignment="left", verticalalignment="top")

    return ax, ax2


def main():
    import matplotlib.pyplot as plt
    
    UPTAKE_RATES_OUTDIR = "parameters/uptake_rates/"
    UPTAKE_RATES_PLOTS_OUTDIR = "out/uptake_rates/"
    DATA_FILE = "data/clean/CUE2/dFBA.pkl"

    # Load condition data
    with open(DATA_FILE, "rb") as f:
        data = pickle.load(f)
    # drawdown_data = pd.read_csv(DRAWDOWN_DATA_CLEAN)

    # Ensure output directories exist
    os.makedirs(UPTAKE_RATES_PLOTS_OUTDIR, exist_ok=True)
    os.makedirs("out/uptake_rates/exponential_regions/", exist_ok=True)
    os.makedirs("out/uptake_rates/Vmaxes/", exist_ok=True)
    os.makedirs("out/uptake_rates/K_M_sanity_check/", exist_ok=True)
    os.makedirs("out/uptake_rates/K_M_fit/", exist_ok=True)

    # Start with Zac's data
    rates = {}
    for (g, a), dat in data.items():
        if g > 0 and a > 0:
            continue

        # Get the data
        t_X_obs = dat["raw"]["raw_b_t"]
        X_obs_replicates = dat["raw"]["raw_b"]
        t_N_obs = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_t"]
        N_obs_replicates = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_s"]

        # Get exponential region
        exponential_region = get_exponential_region(t_X_obs, X_obs_replicates)

        # Plot and save exponential region
        fig, ax = plt.subplots()
        plot_exponential_region_finder(ax, t_X_obs, X_obs_replicates, exponential_region)
        fig.savefig(f"out/uptake_rates/exponential_regions/g{g}_a{a}.png")
        plt.close(fig)

        # Fit Vmax
        fig, ax = plt.subplots()
        mu, N0, Vmax = fit_exp_region_params(t_N_obs, N_obs_replicates, t_X_obs, X_obs_replicates, exponential_region)
        plot_exp_region_params(ax, mu, N0, Vmax, t_N_obs, N_obs_replicates, t_X_obs, X_obs_replicates, exponential_region)
        fig.savefig(f"out/uptake_rates/Vmaxes/g{g}_a{a}.png")
        plt.close(fig)

        # Before fitting K_M, plot sanity check (model error vs. K_M)
        fig, ax = plt.subplots()
        plot_error_vs_K_M(ax, Vmax, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates)
        fig.savefig(f"out/uptake_rates/K_M_sanity_check/g{g}_a{a}.png")
        plt.close(fig)

        # Fit K_M
        K_M = fit_K_M(Vmax, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates)

        # Plot K_M fit
        fig, ax = plt.subplots()
        plot_K_M_fit(ax, Vmax, K_M, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates)
        fig.savefig(f"out/uptake_rates/K_M_fit/g{g}_a{a}.png")
        plt.close(fig)

        # Store fits
        rates["glucose" if g > 0 else "acetate"] = {
            "Vmax": Vmax,
            "K_M": K_M
        }



    # Get data for each metabolite, namely biomass traces (t_X, X) and substrate traces (t_S, S)
    # metabolite_traces = {
    #     metabolite: {
    #         "t_X": growth_data["time (h)"][~growth_data[f"{metabolite}_predicted_mass"].isna()].values * u.hr,
    #         "X": growth_data[f"{metabolite}_predicted_mass"].dropna().values * u.g / COLONY_VOLUME.to("L"),
    #         "t_S": np.array([0, dt_hr]) * u.hr,
    #         "S": np.array([S_0, S_final]) * u.mM
    #     }
    #     for metabolite, dt_hr, S_0, S_final in zip(drawdown_data["Compound"],
    #                                                drawdown_data["dt_hr"],
    #                                                drawdown_data["InitialMetabolite_mM"],
    #                                                drawdown_data["FinalMetabolite_mM"])
    # }

    # Run fitting
    # rates = {
    #     metabolite: fit_M_M_for_given_K_M(trace["t_S"].magnitude,
    #                                       trace["S"].magnitude,
    #                                       trace["t_X"].magnitude,
    #                                       trace["X"].magnitude,
    #                                       K_M.magnitude)
    #     for metabolite, trace
    #     in tqdm(metabolite_traces.items(), desc="Fitting")
    # }

    # # Get raw points for plot
    # substrate_raw_data = {
    #     metabolite: (
    #         raw_data[raw_data["Compound"] == metabolite]["SampleType"].replace(
    #             {"initial": 0, "final": dt_hr}).values * u.hr,
    #         mM_to_peak * raw_data[raw_data["Compound"]
    #                               == metabolite]["IntegratedPeakArea"].values * u.mM
    #     )
    #     for metabolite, dt_hr, mM_to_peak
    #     in tqdm(zip(drawdown_data["Compound"],
    #                 drawdown_data["dt_hr"],
    #                 drawdown_data["mM_to_peak_ratio"]),
    #             desc="Plotting")
    # }
    # substrate_raw_data["acetate"] = (
    #     acetate_data[acetate_data["Type"] ==
    #                  "drawdown (umol)"]["Time (h)"].values * u.hr,
    #     (acetate_data[acetate_data["Type"] ==
    #                   "drawdown (umol)"]["Value"].values * u.umol / CUE_VOLUME).to("mM")
    # )

    # # Make and save plots
    # os.makedirs(UPTAKE_RATES_PLOTS_OUTDIR, exist_ok=True)
    # plot_fits(metabolite_traces, substrate_raw_data,
    #           rates, UPTAKE_RATES_PLOTS_OUTDIR)
    
    # Output rates to json
    os.makedirs(UPTAKE_RATES_OUTDIR, exist_ok=True)
    with open(os.path.join(UPTAKE_RATES_OUTDIR, "fitted_uptake_rates.json"), "w") as f:
        json.dump(rates, f, indent=4)


if __name__ == "__main__":
    main()
