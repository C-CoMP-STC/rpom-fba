import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from cobra.io import read_sbml_model
from tqdm import tqdm
from utils.cobra_utils import get_or_create_exchange, set_active_bound


def main():
    MODEL="model/Rpom_05.xml"
    GROWTH_RATES="parameters/growth_rates/fitted_growth_rates.csv"
    CARBON_SOURCE_IDS="parameters/uptake_rates/carbon_source_ids.json"
    OUTDIR="out/reaction_essentiality"

    os.makedirs(OUTDIR, exist_ok=True)

    model = read_sbml_model(MODEL)

    # Load fitted growth rates
    growth_rates = pd.read_csv(GROWTH_RATES)

    # Load carbon source names -> IDs
    with open(CARBON_SOURCE_IDS, "r") as f:
        carbon_source_ids = json.load(f)

    # Two simultaneous tests:
    # 1) Constrain biomass to produce (at most) the observed growth rate,
    #   see what fluxes are different
    # 2) Remove all reactions one by one, categorize result by
    #   whether it affects biomass and whether it is fatal
    d_flux = np.zeros((len(carbon_source_ids), len(model.reactions)))

    non_boundary_rxns = model.reactions - model.boundary - [model.reactions.get_by_id("RPOM_provisional_biomass")]
    growth_fold_change = np.zeros((len(carbon_source_ids), len(non_boundary_rxns)))
    for r, carbon_source in enumerate(carbon_source_ids.keys()):
        with model:
            carbon_source_id = carbon_source_ids[carbon_source]
            ex = get_or_create_exchange(model, carbon_source_id)
            growth_rate = growth_rates[carbon_source].mean()

            set_active_bound(ex, abs(float(ex.annotation["Experimental rate"])))
            
            baseline_fluxes = model.optimize().fluxes
            biomass = model.reactions.get_by_id("RPOM_provisional_biomass")
            biomass.upper_bound = growth_rate
            constrained_fluxes = model.optimize().fluxes

            d_flux[r, :] = constrained_fluxes - baseline_fluxes
        
        # Experiment 2
        with model:
            # Remove maintenance
            atpm = model.reactions.get_by_id("ATPM")
            atpm.bounds = (0, 0)

            # Set growth on substrate
            carbon_source_id = carbon_source_ids[carbon_source]
            ex = get_or_create_exchange(model, carbon_source_id)
            set_active_bound(ex, abs(float(ex.annotation["Experimental rate"])))

            # Baseline growth
            baseline = model.optimize().objective_value

            # Knockout of every reaction
            for c, rxn in enumerate(tqdm(non_boundary_rxns)):
                # Knock out reaction, get growth, restore
                bounds = rxn.bounds
                rxn.bounds = (0, 0)
                growth_fold_change[r, c] = model.optimize().objective_value
                rxn.bounds = bounds
            
            # Divide by baseline to get fold-change
            if baseline != 0:
                growth_fold_change[r, :] = growth_fold_change[r, :] / baseline
            else:
                growth_fold_change[r, np.nonzero(growth_fold_change[r, :])] = np.inf

    # Exclude zero columns, sort by importance, while keeping associated rxn_ids in order
    rxn_ids = np.array([truncate_str(rxn.id, 15, "...") for rxn in model.reactions])
    nonzero_cols = np.argwhere(~np.all(d_flux[..., :] == 0, axis = 0)).flatten()
    plot_data = d_flux[:, nonzero_cols].T
    rxn_ids = rxn_ids[nonzero_cols]
    
    sorted_order = np.argsort(-np.abs(plot_data).mean(axis=1))
    plot_data = plot_data[sorted_order, :]
    rxn_ids = rxn_ids[sorted_order]
    
    # Plot
    fig, ax = plt.subplots()
    im = ax.imshow(plot_data,
                   cmap="RdYlGn",
                   norm=matplotlib.colors.SymLogNorm(linthresh=0.01),
                   aspect="auto",
                   interpolation='none')
    ax.xaxis.tick_top()
    ax.set_xticks(range(plot_data.shape[1]), labels=carbon_source_ids.keys(), rotation=-90, va="bottom")
    ax.set_yticks(range(plot_data.shape[0]), labels=rxn_ids)
    cb = fig.colorbar(im, ax=ax, extend="both")
    cb.set_label("Less flux under constraint")
    cb.ax.set_title("More flux under constraint")

    fig_size = matplotlib.figure.figaspect(plot_data)
    fig.set_size_inches(2 * fig_size[0], 5 * fig_size[1])
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "constrained_biomass_diffs.png"), dpi=300)

    # Save data to csv
    df = pd.DataFrame({"Reaction" : rxn_ids})
    for substrate, flux_diffs in zip(carbon_source_ids.keys(), plot_data.T):
        df[substrate] = flux_diffs
    df.to_csv(os.path.join(OUTDIR, "constrained_biomass_diffs.csv"), index=False)


    df2 = pd.DataFrame({"Reaction" : [rxn.id for rxn in non_boundary_rxns]})
    for substrate, substrate_data in zip(carbon_source_ids.keys(), growth_fold_change):
        df2[substrate] = substrate_data
    df2.to_csv(os.path.join(OUTDIR, "reaction_essentiality.csv"), index=False)
    

def truncate_str(s, n, suff=None):
    return s if len(s) <= n else (s[:n] + (suff if suff is not None else ""))


if __name__ == "__main__":
    main()
