from pyvis.network import Network
import numpy as np
from cobra.io import read_sbml_model
from experiments.fast_dFBA import setup_drawdown
from tqdm import tqdm


CURRENCY_METABOLITES = {"PROTON[c]",
                        "WATER[c]",
                        "ATP[c]",
                        "NAD[c]",
                        "NADH[c]",
                        "Pi[c]",
                        "CO-A[c]",
                        "ADP[c]",
                        "PPI[c]",
                        "CARBON-DIOXIDE[c]",
                        "NADP[c]",
                        "NADPH[c]"}


def get_graph(model,
              fluxes,
              flux_threshold=0,
              highlight_nodes=[],
              currency=CURRENCY_METABOLITES,
              biomass_color="gray",
              tpp_color="orange"):
    
    # Filter to nonzero fluxes > threshold
    fluxes = fluxes[(fluxes != 0)
                    & (abs(fluxes) >= abs(flux_threshold))].reset_index()
    
    # Normalize fluxes for plotting
    normalize_flux = lambda fx: 5 + 5 * np.log(1 + np.abs(fx))
    fluxes["fluxes"] = normalize_flux(fluxes["fluxes"])

    flux_min = fluxes["fluxes"].min()
    flux_max = fluxes["fluxes"].max()

    g = Network()
    for _, (rxn_id, flux) in tqdm(fluxes.iterrows()):
        # Get reaction object
        rxn = model.reactions.get_by_id(rxn_id)

        # Get color for reaction
        t = (flux - flux_min) / (flux_max - flux_min)
        rxn_color = f"rgba({t * 255}, 0, {(1 - t) * 255}, {0.5 + t/2})"
        if rxn in model.boundary:
            rxn_color = tpp_color
        if rxn_id == "RPOM_provisional_biomass":
            rxn_color = biomass_color
        if rxn_id in highlight_nodes:
            rxn_color = "green"

        # Add reaction node
        g.add_node(rxn_id,
                   label=rxn.build_reaction_string(use_metabolite_names=True),
                   title=rxn_id,
                   value=flux,
                   scaling={"min": flux_min, "max": flux_max},
                   shape="diamond",
                   color=rxn_color)

        # Add nodes for reactants and products
        for met, coeff in rxn.metabolites.items():
            reactant = coeff < 0
            met_id = met.id if met.id not in currency else f"{rxn_id}:{met.id}"

            met_color = "#97c2fc"
            if met.id in currency:
                met_color = "gray"
            if met.id in highlight_nodes:
                met_color = "green"

            g.add_node(met_id,
                       label=f"{met.name}[{met.compartment}]",
                       title=met_id,
                       size=15 if met.id not in currency or met.id in highlight_nodes else 5,
                       color=met_color)
            g.add_edge(met_id if reactant else rxn_id,
                       rxn_id if reactant else met_id,
                       value=flux,
                       scaling={"min": flux_min, "max": flux_max},
                       arrows="to",
                       color=rxn_color)

    return g


def main():
    MODEL = "model/Rpom_05.xml"
    model = read_sbml_model(MODEL)

    with model:
        setup_drawdown(model)
        model.reactions.get_by_id("EX_glc").bounds = (-10, 1000)

        fluxes = model.optimize().fluxes
        g = get_graph(
            model,
            fluxes,
            flux_threshold=fluxes[fluxes > 0].quantile(0.4),
            currency=CURRENCY_METABOLITES - {"CARBON-DIOXIDE[c]"},
            highlight_nodes=["CARBON-DIOXIDE[c]",
                             "1.2.1.2-RXN",
                             "GCVMULTI-RXN",
                             "RXN0-1134",
                             "1.1.1.39-RXN",
                             "RXN-9952",
                             "ISOCITDEH-RXN"])

    g.barnes_hut(gravity=-1000, central_gravity=1)
    g.show_buttons(filter_=['edges', 'physics'])
    g.save_graph("out/g.html")


if __name__ == "__main__":
    main()
