import json
import os
import functools
from abc import abstractmethod
from dataclasses import dataclass
from itertools import pairwise

import cobra
import cobra.util
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra.core.model import Model
from cobra.core.reaction import Reaction
from cobra.io import read_sbml_model
from gem2cue import utils as cue
from tqdm import tqdm

from data.files import DRAWDOWN_DATA_CLEAN, GROWTH_DATA_CLEAN
from parameters.drawdown import *
from parameters.fit_uptake_rates import michaelis_menten_dynamic_system
from utils.cobra_utils import get_or_create_exchange, set_active_bound
from utils.math import euler, get_interpolator, runge_kutta, euler_step, rk_step
from utils.units import u


@dataclass
class Bounds:
    exchange: Reaction
    substrate_id: str

    @abstractmethod
    def bound(self, concentration, biomass, t):
        pass


class ConstantBounds(Bounds):
    def __init__(self, exchange, substrate_id, V, dt=0.01) -> None:
        self.exchange = exchange
        self.substrate_id = substrate_id
        self.V = V
        self.dt = dt

    def bound(self, concentration, biomass, t):
        concentration = max(concentration, 0)
        env_bound = min(0, -concentration / self.dt / biomass)
        bound = max(self.V, env_bound)
        # set_active_bound(self.exchange, bound)
        set_active_bound(self.exchange, bound, abs_bounds=False)


def cache_notify(func):
    func = functools.cache(func)
    def notify_wrapper(*args, **kwargs):
        stats = func.cache_info()
        hits = stats.hits
        results = func(*args, **kwargs)
        stats = func.cache_info()
        if stats.hits > hits:
            print(f"NOTE: {func.__name__}() results were cached")
        return results
    return notify_wrapper


class CachingBoundedOptimizer:
    def __init__(self, model, biomass_id, dynamic_medium):
        self.model = model
        self.biomass_id = biomass_id
        self.dynamic_medium = dynamic_medium
        self.cache = {}
    
    def bound_and_optimize(self,
                           biomass,
                           concentrations,
                           t,
                           use_cache=False):
        with self.model:
            # Bound
            bound_values = []
            for bound, conc in zip(self.dynamic_medium, concentrations):
                bound.bound(conc, biomass, t)
                bound_values.append(bound.exchange.lower_bound)
            
            # Optimize, potentially using cache
            if use_cache:
                key = tuple(bound_values)
                if key in self.cache:
                    # print('hit!')
                    return np.copy(self.cache[key])
            
            # Fall back to calculation
            # print(f"memoizing: {tuple(bound_values)}")
            sol = self.model.optimize()
            exchanges = [bound.exchange for bound in self.dynamic_medium]
            fluxes = np.array([sol.objective_value] + [sol.fluxes[ex.id] for ex in exchanges])
            if use_cache:
                self.cache[key] = np.copy(fluxes)
            
            return np.copy(fluxes)


def bound_and_optimize(model,
                       biomass_id,
                       dynamic_medium,
                       biomass,
                       concentrations,
                       t,
                       callback=None):

    # TODO: Implement caching for bounded_optimize - doesn't work when creating the function each time
    @cache_notify
    def bounded_optimize(exchange_ids, bounds):
        sol = model.optimize()
        fluxes = [sol.objective_value] + [sol.fluxes[ex_id] for ex_id in exchange_ids]
        return fluxes

    exchanges = [bound.exchange for bound in dynamic_medium]

    with model:
        bound(model, biomass_id, dynamic_medium, biomass, concentrations, t)

        # Using lexicographic optimization,
        # first optimize for biomass, then for the exchange fluxes
        # (holding optimal biomass as a constraint),
        # thus guaranteeing a unique optimal set of exchange fluxes.
        
        # lex_constraints = cobra.util.add_lexicographic_constraints(
        #     model,
        #     [biomass_id] + [ex.id for ex in exchanges],
        #     ["max" for _ in range(len(exchanges) + 1)],
        # )
        # fluxes = lex_constraints.values

        # sol = model.optimize()
        # fluxes = [sol.objective_value] + [sol.fluxes[ex.id] for ex in exchanges]
        # print(bound_values)
        fluxes = bounded_optimize(tuple(ex.id for ex in exchanges), tuple(bound_values))

        if callback is not None:
            return fluxes, callback({
                "model": model,
                "biomass_id": biomass_id,
                "dynamic_medium": dynamic_medium,
                "biomass": biomass,
                "concentrations": concentrations,
                "t": t})

    return np.array(fluxes)


def dFBA(
    model,
    biomass_id,
    initial_biomass,
    initial_concentrations,
    dynamic_medium,
    tmax,
    dt=0.01,
    terminate_on_infeasible=True,
    listeners=None,
    desc="",
    integrator="runge_kutta",
    use_cache=False
):
    y0 = np.array([initial_biomass, *initial_concentrations])
    medium_ids = [bounds.exchange.id for bounds in dynamic_medium]

    optimizer = CachingBoundedOptimizer(model, biomass_id, dynamic_medium)

    def df_dt(y, t):
        biomass = y[0]  # g/L
        concentrations = y[1:]  # mM
        
        try:
            fluxes = optimizer.bound_and_optimize(biomass, concentrations, t, use_cache)
            fluxes *= biomass
            return fluxes

        except Exception as e:
            if terminate_on_infeasible:
                raise e
            else:
                return np.array([0 for _ in range(len(medium_ids) + 1)])

    match integrator:
        case "runge_kutta":
            return runge_kutta(
                df_dt,
                y0,
                0,
                tmax,
                dt,
                terminate_on_infeasible,
                pbar_desc=desc,
                listeners=listeners,
            )
        case "euler":
            return euler(
                df_dt,
                y0,
                0,
                tmax,
                dt,
                terminate_on_infeasible,
                pbar_desc=desc,
                listeners=listeners,
            )
        case _:
            raise ValueError(
                f"{integrator} is not a recognized integrator (options: runge_kutta, euler)."
            )


def make_shadow_price_listener(model, substrate_ids, dynamic_medium, n=10):
    def shadow_price_listener(y, t=None):
        substrates = dict(zip(substrate_ids, y[1:]))  # mM

        with model:
            for exchange, bounds in dynamic_medium.items():
                bounds.bound(exchange, substrates[bounds.substrate_id], t=t)

            sol = model.optimize()
            max_shadow_price_metabolites = (
                sol.shadow_prices.abs().sort_values(ascending=False)[:n].index
            )

            return sol.shadow_prices[max_shadow_price_metabolites]

    return shadow_price_listener


def make_cue_listener(model, biomass_id, dynamic_medium, co2_exchange="EX_co2"):
    def cue_listener(y, t=None):
        biomass = y[0]  # * u.g/u.L
        substrates = y[1:]  # mM

        _, (sol, c_ex_rxns) = bound_and_optimize(
            model,
            biomass_id,
            dynamic_medium,
            biomass,
            substrates,
            t,
            callback=lambda d: (model.optimize(), cue.get_c_ex_rxns(d["model"])))

        c_uptake, c_secret = cue.get_c_ex_rxn_fluxes(sol, c_ex_rxns, "cobrapy")
        return cue.calculate_cue(c_uptake, c_secret, co2_exchange)

    return cue_listener


def make_bge_listener(model, biomass_id, dynamic_medium, co2_exchange="EX_co2"):
    def bge_listener(y, t=None):
        biomass = y[0]  # * u.g/u.L
        substrates = y[1:]  # mM

        _, (sol, c_ex_rxns) = bound_and_optimize(
            model,
            biomass_id,
            dynamic_medium,
            biomass,
            substrates,
            t,
            callback=lambda d: (model.optimize(), cue.get_c_ex_rxns(d["model"])))

        c_uptake, c_secret = cue.get_c_ex_rxn_fluxes(sol, c_ex_rxns, "cobrapy")
        return cue.calculate_bge(c_uptake, c_secret, co2_exchange)

    return bge_listener


def make_growth_rate_listener(model, biomass_id, dynamic_medium):
    def growth_rate_listener(y, t=None):
        biomass = y[0]  # * u.g/u.L
        substrates = y[1:]  # mM

        fluxes = bound_and_optimize(
            model,
            biomass_id,
            dynamic_medium,
            biomass,
            substrates,
            t)

        return fluxes[0]

    return growth_rate_listener


def make_boundary_listener(model, biomass_id, dynamic_medium):
    def boundary_listener(y, t=None):
        biomass = y[0]
        substrates = y[1:]

        _, boundary_fluxes = bound_and_optimize(
            model,
            biomass_id,
            dynamic_medium,
            biomass,
            substrates,
            t,
            callback=lambda d: [(rxn, rxn.flux) for rxn in d["model"].boundary if rxn.flux != 0])

        return boundary_fluxes
    
    return boundary_listener


def setup_drawdown(model):
    # Growth is infeasible on the seawater medium as it currently is,
    # needs to be supplemented with FE+2 (also increase everything to 1000 to not be limiting)
    supp_medium = {k: 1000.0 for k, v in model.medium.items()}
    supp_medium["EX_fe2"] = 1000.0
    supp_medium["EX_o2"] = 20.0
    model.medium = supp_medium

    # nadh_dehyd_rxns = [rxn for rxn in model.reactions if rxn.id.startswith("1.6.99.5")]
    # for rxn in nadh_dehyd_rxns:
    #     rxn.add_metabolites({"PROTON[c]" : -4., "PROTON[e]" : 4.})


def plot_data(t, y, carbon_source, initial_C, V_max, t_max, growth_data):
    fig, ax = plt.subplots()

    # Plot data
    ax.plot(
        t,
        (y[:, 0] * u.g / u.L * COLONY_VOLUME.to("L")).to("ug"),
        color="b",
        label="Biomass",
    )
    ax2 = plt.twinx(ax)
    ax2.plot(t, y[:, 1], color="r", label=f"[{carbon_source}] (mM)")

    ax.set_ylabel("Biomass (ug)", color="b")
    ax2.set_ylabel(f"{carbon_source} (mM)", color="r")

    col = f"{carbon_source}_predicted_mass"
    growth_on_carbon_source = growth_data[[col, "time (h)"]]
    growth_on_carbon_source = growth_on_carbon_source[
        ~np.isnan(growth_on_carbon_source[col])
    ]

    ax.plot(
        growth_on_carbon_source["time (h)"],
        (growth_on_carbon_source[col].values * u.g).to("ug"),
        "b--",
        label="Biomass (data)",
    )

    t = growth_data["time (h)"] * u.h
    y = growth_data[f"{carbon_source}_predicted_mass"] * u.g
    mass_curve = get_interpolator(t, y)
    t, x = michaelis_menten_dynamic_system(
        initial_C.magnitude, mass_curve, -abs(V_max), K_M.magnitude, t_max, dt=0.01
    )

    ax2.plot(t, x[:, 0], "r--", label=f"Fitted {carbon_source} drawdown")

    return fig, [ax, ax2]


def plot_shadow_prices(shadow_prices, t):
    shadow_prices = pd.concat(shadow_prices, axis=1,
                              ignore_index=True).T.fillna(0)

    fig, axs = plt.subplots(len(shadow_prices.columns), 1)
    for ax, metabolite in zip(axs, shadow_prices.columns):
        ax.plot(t, shadow_prices[metabolite])
        ax.hlines(0, t.min(), t.max(), colors=["0.4"])

        ax.set_ylabel(metabolite, rotation="horizontal", ha="right")
        ax.xaxis.set_visible(False)
        ax.spines[:].set_visible(False)
        ax.tick_params(left=False, labelleft=False)

    axs[-1].xaxis.set_visible(True)
    axs[-1].spines.bottom.set_visible(True)

    fig.subplots_adjust(hspace=0)
    fig.tight_layout()

    return fig, axs


# def main():
#     MODEL = "model/Rpom_05.xml"
#     BIOMASS_ID = "RPOM_provisional_biomass"
#     CARBON_SOURCES = "parameters/uptake_rates/carbon_source_ids.json"
#     OUTDIR = "out/dFBA/"

#     # Ensure output directory exists
#     os.makedirs(OUTDIR, exist_ok=True)

#     model = read_sbml_model(MODEL)
#     setup_drawdown(model)

#     # Load carbon sources to test
#     with open(CARBON_SOURCES, "r") as f:
#         carbon_sources = json.load(f)
#     carbon_sources = {
#         k: v
#         for k, v in carbon_sources.items()
#         if len(model.metabolites.query(lambda m: m.id == v)) == 1
#     }

#     # Load data for compoarison
#     uptake_data = pd.read_csv(DRAWDOWN_DATA_CLEAN)
#     growth_data = pd.read_csv(GROWTH_DATA_CLEAN)

#     for carbon_source, carbon_source_id in carbon_sources.items():
#         with model:
#             # Get Metabolite object and exchange reaction for the given carbon source
#             exchange_rxn = get_or_create_exchange(model, carbon_source_id)
#             V_max = abs(float(exchange_rxn._annotation["Experimental rate"]))

#             # Get initial conditions
#             initial_C = uptake_data[uptake_data["Compound"] ==
#                                     carbon_source]["InitialMetabolite_mM"].values[0] * u.mM
#             initial_biomass = growth_data[f"{carbon_source}_predicted_mass"].values[0] * u.g / (
#                 COLONY_VOLUME.to("L"))
#             y0 = np.array([initial_biomass.magnitude,
#                            initial_C.magnitude])

#             # Run experiment
#             tmax = uptake_data[uptake_data["Compound"] ==
#                                carbon_source]["dt_hr"].values[0]

#             dynamic_medium = {exchange_rxn: MichaelisMentenBounds(
#                 carbon_source_id, V_max, K_M)}
#             t, y, listeners = dFBA(model, BIOMASS_ID, [carbon_source_id], dynamic_medium, COLONY_VOLUME,
#                                    y0, tmax, dt=0.01, terminate_on_infeasible=True,
#                                    listeners=[make_shadow_price_listener(
#                                        model, [carbon_source_id], dynamic_medium)],
#                                    desc=carbon_source)

#             # Plot shadow prices over time
#             fig_sp, axs_sp = plot_shadow_prices(listeners, t)
#             axs_sp[0].set_title(carbon_source)
#             fig_sp.savefig(os.path.join(
#                 OUTDIR, f"{carbon_source} shadow prices.png"), bbox_inches='tight')

#             # Plot data
#             fig, _ = plot_data(t, y, carbon_source, initial_C,
#                                V_max, tmax, growth_data)
#             fig.set_size_inches(5, 3)
#             fig.tight_layout()
#             fig.savefig(os.path.join(OUTDIR, f"{carbon_source} dFBA.png"))


# if __name__ == "__main__":
#     main()
