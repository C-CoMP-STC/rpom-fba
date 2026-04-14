import marimo

__generated_with = "0.23.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    from cobra.io import load_model

    return load_model, pd


@app.cell
def _(load_model, pd):
    # Load model
    ecoli = load_model("iJO1366")

    # Load kcats, subset to E. coli wild-type
    kcats = pd.read_csv("notebooks/data/gecko_kcats_table.csv")
    kcats = kcats[kcats["record_organism"] == "Escherichia coli"]
    kcats = kcats[
        kcats["record_enzyme_type"].str.contains("wild").values.astype(bool)
        & ~kcats["record_enzyme_type"].str.contains("mutant").values.astype(bool)
    ]

    # Load mapping of ecocyc gene ids to accessions, reactions
    ecocyc_genes_to_accessions = pd.read_csv("notebooks/data/EcoCyc_genes_to_reactions.tsv", sep="\t")
    ecocyc_genes_to_accessions["Reactions of gene"] = ecocyc_genes_to_accessions["Reactions of gene"].str.split(" // ")

    # Load Schmidt proteome dataset, merge in accessions
    prot = pd.read_csv("notebooks/data/schmidt2015_javier_table.tsv", sep="\t")
    prot = prot.merge(ecocyc_genes_to_accessions, left_on="EcoCycID", right_on="Gene Name")
    prot
    return ecoli, kcats, prot


@app.cell
def _(kcats):
    kcats
    return


@app.cell
def _(ecoli, prot):
    rxn_to_enzyme_counts = []

    for rxn in ecoli.reactions:
        matches = prot[prot["Accession-1"].isin([str(g) for g in rxn.genes])]


    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Introduction

    ## Preliminaries

    Metabolic fluxes are subject to a fundamental constraint in the maximal turnover numbers of the associated enzymes.

    Let $v \in \mathbb R_{\ge 0}^{m}$ be the vector of (nonnegative) reaction fluxes, and $K \in \mathbb R ^{m \times p}$ be the flux-enzyme mapping matrix, such that $K_{ij} > 0$ iff enzyme $j$ catalyzes reaction $i$. Nonzero entries of $K$ represent the turnover numbers ($k_{cat}$'s) associated with the corresponding enzyme-reaction pair. Thus we have:

    $$
    v \le K [E],
    $$

    where $[E] \in \mathbb R^p$ is the vector of enzyme concentrations (units?). We can make this inequality exact by introducing a vector of enzyme saturations $f \in \mathbb [0, 1]^{p}$,

    $$
    v = Kf\odot [E]
    $$

    (where $\odot$ is element-wise multiplication, equivalent to $v=K\text{diag}(f)[E]$). We can think of $f$ as representing what fraction of enzymes are "occupied" at a given time, which is often given by Michaelis-Menten-like expressions.

    ## Problem Statement

    Entries of $K$ are difficult to measure, and often vary widely between measurements. Entries of $f$ depend on metabolite concentrations, which are also relatively difficult to measure, and may have unknown functional relationships with the actual values of $f$ (though there is evidence that they often are close to $1$ in certain settings).

    On the other hand, $[E]$ directly comes from proteomic measurements, which exist for many different organisms in varying environmental conditions.

    Given the uncertainty in these key parameters, can we use a Bayesian approach to estimate or regularize measurements of $k_{cats}$ for metabolic enzymes, from measured proteomes?

    ## Model

    Let us introduce the key constraint that under steady-state growth, $S v = 0$, where $S$ is the stoichiometric matrix defined exactly as in FBA. Then we should have

    $$
    S K f \odot [E] = 0
    $$

    in steady-state growth. But practically speaking, under experimental noise we'll instead have

    $$
    S K f \odot [E] \sim \varepsilon
    $$

    where $\varepsilon$ is a zero-mean noise vector (NEED TO COME BACK AND FIGURE OUT WHAT DISTRIBUTION TO USE, CONSIDERING LOGNORMAL K and GAUSSIAN(?) [E]).

    Here, $S$ is known from a genome-scale metabolic model (GEM), and $[E]$ comes from a measured proteome. This leaves $K$ and $f$. We can put prior distributions on these, possibly informed by measurements, and sample their posteriors

    $$
    K_{ij} \sim \text{LogNormal}(\mu_{ij}, \sigma^2_{ij}) \qquad \forall i, j\text{ with measured values}\\
    K_{ij} \sim \text{LogNormal}(\mu_{marg}, \sigma^2_{marg}) \qquad \forall i, j\text{ with known enzyme-reaction relationship}\\
    f_i \sim \text{Beta} (\alpha, \beta) \qquad \forall i
    $$
    """)
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Plan:

    1. Check how well $f=1$ captures _E. coli_ dataset
        1. Take S
    """)
    return


if __name__ == "__main__":
    app.run()
