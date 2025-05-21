import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from utils.units import u
from utils.colors import *

from parameters.drawdown import (
    MASS_PER_CELL,
    MASS_PER_CELL_GLUCOSE,
    MASS_PER_CELL_ACETATE,
    CUE_VOLUME)

C_PER_GLUCOSE = 6
C_PER_ACETATE = 2
DATA_FILE = "data/clean/CUE2/cue2.xlsx"
# Save data for dFBA
OUTFILE_DFBA = "data/clean/CUE2/dFBA.pkl"


def plot_data(ax, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates, lw=2, N_color = C_GLUCOSE):
    ax.plot(t_X_obs, X_obs_replicates.T, marker="o", lw=0, color=C_BIOMASS, alpha=0.25)
    ax.plot(t_X_obs, X_obs_replicates.mean(axis=0), color=C_BIOMASS)
    ax2 = ax.twinx()
    ax2.plot(t_N_obs, N_obs_replicates.T, marker="o", lw=lw, color=N_color)

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Biomass (g/L)")
    ax2.set_ylabel("Substrate (mM)")

    return ax, ax2


def extract_growth_parameters(data, NO_GROWTH_THRESHOLD=0.01,
                              CV_THRESHOLD=0.35, DETECTION_LIMIT=1e-4):
    
    # Create plot output directory
    os.makedirs("data/clean/CUE2/logs", exist_ok=True)

    exponential_regions = {}
    growth_rates = {}
    for (g, a), dat in data.items():
        if g > 0 and a > 0:
            continue

        # Get raw data, putting biomass on log-scale
        t_X_obs = dat["raw"]["raw_b_t"]
        X_obs_replicates = np.log(dat["raw"]["raw_b"] + DETECTION_LIMIT)

        # Plot dlog(X)/dt
        fig, ax = plt.subplots()
        mean_dX_dt = np.mean([np.gradient(X_obs, t_X_obs) for X_obs in X_obs_replicates], axis=0)
        for X_obs in X_obs_replicates:
            dX_dt = np.gradient(X_obs, t_X_obs)
            ax.hlines(0, 0, t_X_obs[-1], color="0.5", linestyle="-", zorder=-1)
            ax.plot(t_X_obs, dX_dt, marker="o", lw=0, color=C_BIOMASS, alpha=0.25)
        ax.plot(t_X_obs, mean_dX_dt, color=C_BIOMASS)
        ax.set_title(f"{g.magnitude:.2f} glc, {a.magnitude:.2f} ace")
        ax.set_xlim(0, t_X_obs[-1])
        ax.set_ylim(-0.1, 0.8)
        ax.set_ylabel(r"$\frac{{d}}{{dt}}\log (X)$", color=C_BIOMASS)
        ax.set_xlabel("Time (h)")

        # Get the standard deviation of the gradient
        ax2 = ax.twinx()
        dX_dt_std = np.std([np.gradient(X_obs, t_X_obs) for X_obs in X_obs_replicates], axis=0)
        dX_dt_CV = dX_dt_std / np.abs(mean_dX_dt)
        ax2.plot(t_X_obs, dX_dt_CV, color="r")
        ax2.set_ylim(0, 1)
        ax2.set_ylabel(r"$\sigma/\mu$", color="r")

        # Okay, let's say that the exponential region is where the standard deviation
        # is less than 0.05, and the gradient is greater than 0.1
        # exponential_region = (dX_dt_std < STD_THRESHOLD) & (mean_dX_dt > NO_GROWTH_THRESHOLD)
        exponential_region = (dX_dt_CV < CV_THRESHOLD) & (mean_dX_dt > NO_GROWTH_THRESHOLD)
        t_exponential = t_X_obs[exponential_region]

        # Store the exponential region
        exponential_regions[(g, a)] = exponential_region

        # Calculate the growth rate
        growth_rate = np.mean(mean_dX_dt[exponential_region])
        growth_rates[(g, a)] = growth_rate

        # Plot the exponential region
        ax.fill_between(t_X_obs, -1, 1, where=exponential_region, color="0.9", zorder=-2)
        ax.hlines(growth_rate, t_exponential[0], t_exponential[-1], color="k", linestyle="--", zorder=-1)
        ax.text(t_exponential[0] + 2, growth_rate + 0.1, f"μ = {growth_rate:.2f}", horizontalalignment="left")

        fig.tight_layout()
        fig.set_size_inches(4, 3)

        fig.savefig(f"data/clean/CUE2/logs/exponential_region_{g.magnitude:.2f}_{a.magnitude:.2f}.png", dpi=300)
        plt.close(fig)

    return exponential_regions, growth_rates


def main():
    cells = pd.read_excel(DATA_FILE, "Cell number", skiprows=1)
    glucose = pd.read_excel(DATA_FILE, "Glucose", skiprows=1)
    acetate = pd.read_excel(DATA_FILE, "Acetate", skiprows=1)
    biomass_carbon = pd.read_excel(DATA_FILE, "C content", skiprows=1)
    carbon = pd.read_excel(DATA_FILE, "CO2 cumulative", skiprows=1)
    bge = pd.read_excel(DATA_FILE, "BGE", skiprows=1)


    conditions = {
        "Glucose": (2.0 * u.mM, 0.0 * u.mM),
        "Acetate": (0.0 * u.mM, 6.0 * u.mM),
        "4gluc_8ace": ((4/C_PER_GLUCOSE) * u.mM, (8/C_PER_ACETATE) * u.mM),
        "8gluc_4ace": ((8/C_PER_GLUCOSE) * u.mM, (4/C_PER_ACETATE) * u.mM),
    }

    data = {}
    for pref, key in conditions.items():
        g, a = key
        # Estimate the mass per cell by interpolating between glucose and acetate measurements
        # This is a very rough estimate, but it's the best we can do with the data we have
        theta = (a / (g + a)).magnitude
        mass_per_cell = np.interp(theta, [0, 1], [MASS_PER_CELL_GLUCOSE.magnitude, MASS_PER_CELL_ACETATE.magnitude]) * u.pg

        # Biomass
        cell_data = cells[[col for col in cells.columns if col.startswith(pref)]]  # Cells / mL
        b_t = cells["Time"].values
        cell_raw = cell_data[[
            col for col in cell_data.columns if not col.endswith("mean")]]
        cell_mean = cell_data[[
            col for col in cell_data.columns if col.endswith("mean")]]
        # Convert to g/L
        raw_b = (cell_raw.values * 1/u.mL * mass_per_cell).to("g/L").magnitude.T
        b_s = (cell_mean.values.T[0] * 1/u.mL * mass_per_cell).to("g/L").magnitude

        # Glucose
        gluc_data = glucose[[
            col for col in glucose.columns if col.startswith(pref)]]
        g_t = glucose["Time (single)" if pref in {
            "Glucose", "Acetate"} else "Time (double)"]
        g_t = g_t[~np.isnan(g_t)]
        glucose_raw = gluc_data.values[:g_t.size, :].T
        glucose_mean = glucose_raw.mean(axis=0)

        # Acetate
        ace_data = acetate[[
            col for col in acetate.columns if col.startswith(pref)]]
        a_t = acetate["Time (single)" if pref in {
            "Glucose", "Acetate"} else "Time (double)"]
        a_t = a_t[~np.isnan(a_t)]
        acetate_raw = ace_data.values[:a_t.size, :].T
        acetate_mean = acetate_raw.mean(axis=0)

        # BGE
        bge_data = bge[[col
                        for col in bge.columns
                        if col.startswith(pref)
                        and not col.endswith("mean")]]
        bge_t = bge["OD time"].values
        bge_raw = bge_data.values.T

        # Remove NaNs
        # nan_pos = np.isnan(bge_raw).any(axis=0)
        # bge_raw = bge_raw[:, ~nan_pos]
        # bge_t = bge_t[~nan_pos]

        bge_mean = bge_raw.mean(axis=0)


        data[key] = {
            "raw": {
                "raw_b_t": b_t,
                "raw_b": raw_b,
                "raw_g_t": g_t,
                "raw_g_s": glucose_raw,
                "raw_a_t": a_t,
                "raw_a_s": acetate_raw,
                "raw_bge_t": bge_t,
                "raw_bge": bge_raw,
            },
            "mean": {
                "g_t": g_t,
                "g_s": glucose_mean,
                "a_t": a_t,
                "a_s": acetate_mean,
                "b_t": b_t,
                "b_s": b_s,
                "bge_t": bge_t,
                "bge": bge_mean
            }
        }

    # Run time-rescaling routine
    # (basically, nutrient and growth curve data were not paired, and
    # a temperature difference resulted in slower growth for the nutrient data.
    # To align the data, we need to rescale the time for the nutrient data)
    exponential_regions, growth_rates = extract_growth_parameters(data)

    alphas = {}
    for (g, a), dat in data.items():
        if g > 0 and a > 0:
            continue

        # Get known variables
        N0 = max(g, a).magnitude
        t_X_obs = dat["raw"]["raw_b_t"]
        t_lag = t_X_obs[exponential_regions[(g, a)]].min()
        t_stationary = t_X_obs[exponential_regions[(g, a)]].max()
        X0 = dat["mean"]["b_s"][dat["mean"]["b_t"] == t_lag][0]
        mu = growth_rates[(g, a)]


        # Get timepoints of N observations
        t_N_obs = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_t"].values

        # Define function to get N(t) values at observation times
        # for a given k, alpha
        def N_t(k, alpha, t=t_N_obs):
            after_t_l =  N0 + ((k * X0) / (mu * alpha)) * (1 - np.exp(mu * (alpha * t - t_lag)))
            result = np.array([N0] * len(t))
            result[t > t_lag / alpha] = after_t_l[t > t_lag / alpha]
            return result
        
        # Optimize k and alpha
        def error(params, lamb=0.5):
            # Term 1: SSE of N(t) vs. N_obs UP UNTIL STATIONARY
            k, alpha = params
            N_t_pred = N_t(k, alpha)
            N_t_pred[N_t_pred < 0] = 0

            N_obs_replicates = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_s"]

            term1 = np.sum((N_obs_replicates - N_t_pred) ** 2)

            # Term 2: penalty for zero point != t_stationary / alpha
            t_zero = t_lag / alpha + (1 / (alpha * mu)) * np.log(1 + N0 * mu * alpha / (k * X0))
            term2 = (t_zero - t_stationary / alpha) ** 2
            return term1 + lamb * term2

        # Plot error surface
        k_vals = np.linspace(1e-6, 20, 100)
        alpha_vals = np.linspace(1e-6, 2, 100)
        K, A = np.meshgrid(k_vals, alpha_vals)
        error_surface = np.zeros(K.shape)
        for i in range(K.shape[0]):
            for j in range(K.shape[1]):
                k = K[i, j]
                alpha = A[i, j]
                error_surface[i, j] = error([k, alpha])
        error_surface = np.log(error_surface)
        fig, ax = plt.subplots()
        ax.contourf(K, A, error_surface, levels=50, cmap="viridis")
        ax.set_xlabel("k")
        ax.set_ylabel("alpha")
        ax.set_title(f"{g.magnitude:.2f} glc, {a.magnitude:.2f} ace")
        fig.colorbar(ax.contourf(K, A, error_surface, levels=50, cmap="viridis"), label="log(error)")
        fig.tight_layout()
        fig.set_size_inches(4, 3)
        fig.savefig(f"data/clean/CUE2/logs/error_surface_{g.magnitude:.2f}_{a.magnitude:.2f}.png", dpi=300)

        minimizer = minimize(error, [1, 1], bounds=[(0, None), (1e-6, None)])
        k, alpha = minimizer.x

        # Store the parameters
        alphas[(g, a)] = alpha

        print(f"Condition: {g.magnitude:.2f} glc, {a.magnitude:.2f} ace")
        print(f"  k: {k:.2f}, alpha: {alpha:.2f}")

        # Plot the data
        fig, ax = plt.subplots()
        t_X_obs = dat["raw"]["raw_b_t"]
        X_obs_replicates = dat["raw"]["raw_b"]
        t_N_obs = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_t"]
        N_obs_replicates = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_s"]
        ax, ax2 = plot_data(ax, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates, lw=0.5, N_color = C_GLUCOSE if g > 0 else C_ACETATE)

        # Shade exponential region
        ax.fill_between(t_X_obs, 0, X_obs_replicates.max(), where=exponential_regions[(g, a)], color="0.9", zorder=-2)

        # Plot the fitted N(t)
        N_t_pred = N_t(k, alpha, t=t_X_obs)
        N_t_pred[N_t_pred < 0] = 0
        ax2.plot(t_X_obs, N_t_pred, color="k", lw=2)

        ax.set_title(f"{g.magnitude:.2f} glc, {a.magnitude:.2f} ace")

        fig.tight_layout()
        fig.set_size_inches(4, 3)
        fig.savefig(f"data/clean/CUE2/logs/N_t_fit_{g.magnitude:.2f}_{a.magnitude:.2f}.png", dpi=300)
        plt.close(fig)

        # Store the time-rescaled data
        t_N_obs = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_t"]
        t_N_obs_rescaled = t_N_obs * alpha
        data[(g, a)]["raw"][f"raw_{'g' if g > 0 else 'a'}_t"] = t_N_obs_rescaled

        # Plot after refitting
        fig, ax = plt.subplots()
        t_X_obs = dat["raw"]["raw_b_t"]
        X_obs_replicates = dat["raw"]["raw_b"]
        t_N_obs = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_t"]
        N_obs_replicates = dat["raw"][f"raw_{'g' if g > 0 else 'a'}_s"]
        ax, ax2 = plot_data(ax, t_X_obs, X_obs_replicates, t_N_obs, N_obs_replicates, lw=0.5, N_color = C_GLUCOSE if g > 0 else C_ACETATE)
        # Shade exponential region
        ax.fill_between(t_X_obs, 0, X_obs_replicates.max(), where=exponential_regions[(g, a)], color="0.9", zorder=-2)
        ax.set_title(f"{g.magnitude:.2f} glc, {a.magnitude:.2f} ace")
        fig.tight_layout()
        fig.set_size_inches(4, 3)
        fig.savefig(f"data/clean/CUE2/logs/N_t_fit_rescaled_{g.magnitude:.2f}_{a.magnitude:.2f}.png", dpi=300)
        plt.close(fig)


    # Save data for dFBA
    with open(OUTFILE_DFBA, "wb") as f:
        pickle.dump(data, f)

if __name__ == "__main__":
    main()