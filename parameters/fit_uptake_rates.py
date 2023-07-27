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
METHOD_NOT_RECOGNIZED = "{} is not a recognized method to fit uptake rates!"
BG_GRAY = "0.4"

# Dry mass
INITIAL_MASS = 405.3053598e-6  # g

# Michaelis-Menten
# Assuming K_M fixed for all metabolites, consider varying later
K_M = 1


def get_mass_interpolator(metabolite, growth_data):
    mass = growth_data[f"{metabolite}_predicted_mass"]
    mass = mass[~np.isnan(mass)].values
    t_range = growth_data["time (h)"][:len(mass)].values

    def X_t(t): return np.interp([t], t_range, mass)[0]

    return X_t


def michaelis_menten_dynamic_system(S_0, X_t, V_max, K_M, max_t, dt=0.1):
    time = np.arange(0, max_t, dt)
    trajectory = np.zeros((time.size, 2))

    s = S_0
    x = X_t(0)
    for i, t in enumerate(time):
        trajectory[i, :] = [s, x]

        x = X_t(t)
        dS = V_max * (s / (K_M + s)) * x

        s += dS * dt

    return np.array(time), np.array(trajectory)


def fit_M_M_for_given_K_M(S_0, X_t, S_final, K_M, max_t, dt=0.1):
    result = minimize(
        (
            lambda v_max: (
                S_final -
                michaelis_menten_dynamic_system(
                    S_0, X_t, v_max, K_M, max_t, dt=dt)[1][-1, 0]
            )**2
        ),
        [1])
    return result


def fit_uptake_rates(
    data,
    mass_interpolators,
    growth_rates="GrowthRate",
    metabolites="Compound",
    initial_conc="InitialMetabolite_mM",
    final_conc="FinalMetabolite_mM",
    dts="dt_hr",
    method="michaelis-menten"
):
    rates = {}

    # Prevent taking log of 0 using this method:
    def correct_zeros(x): return x if x != 0 else .00001

    match method:
        case "michaelis-menten":
            rates = {
                metabolite: fit_M_M_for_given_K_M(initial_met,
                                                  mass_interpolators[metabolite],
                                                  final_met,
                                                  K_M,
                                                  tmax).x[0]
                for metabolite, final_met, initial_met, tmax
                in zip(data[metabolites],
                       data[final_conc],
                       data[initial_conc],
                       data[dts])
            }
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
    mass_interpolators,
    raw_data,
    rates,
    metabolites="Compound",
    dt="dt_hr",
    growth_rate="GrowthRate",
    method="michaelis-menten"
):
    # Ensure output directory exists
    os.makedirs(UPTAKE_RATES_PLOTS_OUTDIR, exist_ok=True)

    for metabolite, rate in rates.items():
        # Subset data to just one metabolite
        metabolite_data = raw_data[raw_data[metabolites] == metabolite]

        total_time = data[data[metabolites] == metabolite][dt].values[0]
        mu = data[data[metabolites] == metabolite][growth_rate].values[0]
        mM_to_peak_ratio = data[data[metabolites] ==
                                metabolite]["mM_to_peak_ratio"].values[0]

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
            case "michaelis-menten":
                t, x = michaelis_menten_dynamic_system(
                    t0.mean(), mass_interpolators[metabolite], rate, K_M, total_time)

                ax.plot(t, x[:, 0], color="g", alpha=0.5, label=f"[S]")
                ax.set_xlabel("Time (hr)")
                ax.set_ylabel("Substrate concentration (mM)")
                ax2 = ax.twinx()
                ax2.plot(t, x[:, 1], color="r", label="X(t)")
                ax2.set_ylabel("Colony mass (g)", color="r")
            case "linear":
                ax.plot(t,
                        mean_initial + (rate * INITIAL_MASS * (np.exp(mu * t) - 1) / mu))
            case "exponential":
                ax.plot(t, mean_initial * np.exp(rate * t))
            case _:
                raise ValueError(METHOD_NOT_RECOGNIZED.format(method))

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

    growth_data = pd.read_csv(GROWTH_DATA_FILE)
    mass_interpolators = {
        metabolite: get_mass_interpolator(metabolite, growth_data)
        for metabolite in data["Compound"]
    }

    growth_rates = (pd.read_csv(GROWTH_RATES_FILE)
                    .mean()
                    .reset_index(name="GrowthRate")
                    .rename(columns={"index": "Compound"}))
    data = data.merge(growth_rates, how="inner", on="Compound", validate="1:1")
    rates = fit_uptake_rates(data, mass_interpolators,
                             method="michaelis-menten")

    print("Plotting fitted rates...")
    raw_data = pd.read_excel(RAW_UPTAKE_RATES_FILE,
                             sheet_name=RAW_UPTAKE_RATES_SHEET)
    plot_fits(data, mass_interpolators, raw_data, rates)

    # Output rates to json
    print("Saving fitted rates to json...")
    os.makedirs(UPTAKE_RATES_OUTDIR, exist_ok=True)

    # Writing to .json
    outfile = os.path.join(UPTAKE_RATES_OUTDIR, "fitted_uptake_rates.json")
    with open(outfile, "w") as f:
        json.dump(rates, f, indent=4)


if __name__ == "__main__":
    main()
