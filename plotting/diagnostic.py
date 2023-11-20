import numpy as np
import matplotlib.pyplot as plt


def formula_color(elements):
    element_colors = {
        "C" : np.array([0, 0, 1]),
        "N" : np.array([0, 1, 0]),
        "P" : np.array([1, 0, 1]),
        "H" : np.array([1, 0, 0]),
        "O" : np.array([0, 1, 0.5])
    }

    denom = sum(elements.values())

    result = sum(
        (elements.get(elem, 0) / denom) * c
        for elem, c in element_colors.items()
    )
    if any(result > 1):
        result /= result.sum()
    return result


def uptakes_and_secretions(model, boundary_fluxes, sol, initial_biomass, CUE_VOLUME, initial_glucose, initial_acetate):
    fig, ax = plt.subplots()

    sum_secreting = 0
    sum_uptake = 0
    for _, (rxn, flux, secreting, molmass) in boundary_fluxes.iterrows():
        flux = abs(flux) * molmass  # flux is in mmol/hr, so multiplying by molmass results in mg/hr
        if flux == 0:
            continue

        ax.bar(int(secreting),
            flux,
            width=0.8,
            bottom=sum_secreting if secreting else sum_uptake,
            align="center",
            color=formula_color(list(model.reactions.get_by_id(rxn).metabolites.keys())[0].elements),
            label=str(rxn))

        if secreting:
            sum_secreting += flux
        else:
            sum_uptake += flux

    ax.bar(1,
            sol.objective_value * initial_biomass * CUE_VOLUME.to("L").magnitude * 1000,
            width=0.8,
            bottom=sum_secreting,
            color="gray",
            align="center",
            label="biomass")
    # ax.hlines(0, 0, 1, ["k"])
    ax.set_xticks([0, 1], ["uptake", "secretions"])
    ax.set_ylabel("flux (mg / hr)")
    ax.set_title(f"{initial_glucose} mM Glucose, {initial_acetate} mM Acetate")
    fig.legend()

    return fig, ax



def main():
    pass

if __name__ == "__main__":
    main()
