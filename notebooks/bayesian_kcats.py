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
    import pymc as pm
    import arviz as az

    from cobra.io import load_model
    from scipy.stats import norm, lognorm, beta

    return az, beta, load_model, lognorm, norm, np, pd, plt, pm


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
    return ecocyc_genes_to_accessions, ecoli, kcats, prot


@app.cell
def _(kcats):
    kcats
    return


@app.cell
def _(ecocyc_genes_to_accessions):
    ecocyc_genes_to_accessions
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


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Non-probabilistic analysis

    In the simplest case, let's first $f=\mathbf 1$ and assume some $K_{ij}$'s are measured:

    $$
    S \underbrace{(\hat K + K^?)}_{K} [E] \approx 0
    $$

    where we've decomposed $K$ into $\hat K$, the matrix of measured $k_{cat}$'s, and $K^?$, the unknown $k_{cat}$ values. Both these are $\mathbb R_{\ge 0}^{m \times p}$ matrices, and while we don't know the entries of $K^?$, we do know which entries are non-zero.

    In particular, some columns
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


@app.cell
def _(ecoli):
    ecoli.genes
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Toy PyMC example

    Let's take the following simple bifurcating and merging reaction network:

    $$
    \begin{align*}
    \emptyset \stackrel{v_1}\rightarrow A\\
    A \stackrel{v_2}\rightarrow B\\
    A \stackrel{v_3}\rightarrow C\\
    B \stackrel{v_4}\rightarrow X\\
    C \stackrel{v_5}\rightarrow X\\
    X \stackrel{v_6}\rightarrow \emptyset.
    \end{align*}
    $$

    We'll fix $v_1=10$, and observe $k_{cat}$'s for $v_4$ and $v_5$, leaving $k_{cat}$'s for $v_1, v_2, v_3,$ and $v_6$ unknown. We'll observe protein concentrations $E_i$ for all reactions. Ideally, as we vary $k_{cat}^{(4)}$ and $k_{cat}^{(5)}$, we should see $k_{cat}^{(2)}$ and $k_{cat}^{(3)}$ vary accordingly.

    Our stoichiometric matrix and bounds vectors are, respectively:

    $$
    \begin{align*}
    S &= \begin{pmatrix}
        1 & -1 & -1 &    &    &   \\
          &  1 &    & -1 &    &   \\
          &    &  1 &    & -1 &   \\
          &    &    &  1 &  1 & -1
    \end{pmatrix},\\
    \ell &= \begin{pmatrix} 10 & 0 & 0 & 0 & 0 & 0\end{pmatrix}\\
    u &=    \begin{pmatrix} 10 & - & - & - & - & -\end{pmatrix}
    \end{align*}
    $$

    (where the omitted values in $u$ represent some suitably large upper bounds).

    Each reaction will further have two associated random variables, $k_{cat}^{(i)}$ and $f_i$. For unobserved $k_{cat}$'s, we'll start with the uninformative prior $\log k_{cat}^{(i)} \sim \mathcal N(0.5, 2.5^2)$, informed by the marginal $k_{cat}$ distribution from above. For the $f_i$'s, we'll use $\text{Beta}$ priors. We'll first consider the case of a single shared $f$, then see how identifiability may be lost with individual $f_i$'s per reaction. Observed $k_{cat}$'s will be treated as deterministic.
    """)
    return


@app.cell(hide_code=True)
def _(kcats, lognorm, norm, np, plt):
    # Fit marginal distribution for log-kcats
    logkcat_marginal = norm(*norm.fit(np.log(kcats["effective_kcat_s"])))
    kcat_marginal = lognorm(s = logkcat_marginal.std(), scale = np.exp(logkcat_marginal.mean()))

    _fig, (_ax1, _ax2) = plt.subplots(1, 2)

    _ax1.hist(np.log(kcats["effective_kcat_s"]), bins=50, density=True)
    _x = np.linspace(np.log(kcats["effective_kcat_s"]).min(), np.log(kcats["effective_kcat_s"]).max(), 100)
    _ax1.plot(_x, logkcat_marginal.pdf(_x))
    _ax1.set_xlabel(r"$\log k_{cat}$")
    _ax1.set_ylabel("Density")
    _ax1.set_title(rf"Marginal $\log k_{{cat}} \sim N ({logkcat_marginal.mean():.2f}, {logkcat_marginal.std():.2f})$")

    _ax2.hist(kcats["effective_kcat_s"], bins=50, density=True)
    _x = np.linspace(kcats["effective_kcat_s"].min(), kcats["effective_kcat_s"].max(), 100)
    _ax2.plot(_x, kcat_marginal.pdf(_x))
    _ax2.set_xlabel(r"$k_{cat}$")
    _ax2.set_yscale("log")
    _ax2.set_title(r"Marginal $k_{{cat}} \sim \log N(\cdot)$")

    _fig.set_size_inches(8,2.5)
    _fig.tight_layout()
    _fig
    return kcat_marginal, logkcat_marginal


@app.cell(hide_code=True)
def _(mo):
    alpha_1 = mo.ui.slider(0.01, 10, 0.01, value=5, label=r"$\alpha$", show_value=True)
    beta_1 = mo.ui.slider(0.01, 10, 0.01, value=1, label=r"$\beta$", show_value=True)
    return alpha_1, beta_1


@app.cell(hide_code=True)
def _(alpha_1, beta, beta_1, mo, np, plt):
    f_prior_1 = beta(a=alpha_1.value, b=beta_1.value)

    _fig, _ax = plt.subplots()

    _x = np.linspace(0, 1, 100)
    _y = f_prior_1.pdf(_x)
    _ax.plot(_x, _y)
    _ax.set_ylim(-0.1, _y[~np.isnan(_y) & np.isfinite(_y)].max())
    _ax.set_xlabel("$f$")
    _ax.set_ylabel("Density")

    _fig.set_size_inches(3, 1.5)
    _fig.tight_layout()

    mo.vstack([alpha_1, beta_1, _fig], heights="equal")
    return (f_prior_1,)


@app.cell(hide_code=True)
def _(alpha_1, beta_1, logkcat_marginal, np, pm):
    # Stoichiometric matrix
    _S = np.array([
        [1, -1, -1, 0, 0, 0],
        [0, 1, 0, -1, 0, 0],
        [0, 0, 1, 0, -1, 0],
        [0, 0, 0, 1, 1, -1]
    ])

    # Observed enzyme concentrations
    _E = np.array([1, 1, 1, 1, 1, 1])

    # Observed kcats
    _observed_kcats = np.array([10, 20])
    _observed_kcat_indices = np.array([3, 4])

    # Observed (fixed) fluxes
    _observed_fluxes = np.array([10])
    _observed_flux_indices = np.array([0])

    _coords = {
        "reactions": np.arange(_S.shape[1]),
        "metabolites": np.arange(_S.shape[0]),
    }
    with pm.Model(coords=_coords) as model:
        # Priors for log-kcats
        _logkcat = pm.Normal("log_kcat", mu=logkcat_marginal.mean(), sigma=logkcat_marginal.std(), dims="reactions")
        _kcat = pm.Deterministic("kcat", np.exp(_logkcat), dims="reactions")

        # Prior for saturation f
        _f = pm.Beta("f", alpha=alpha_1.value, beta=beta_1.value)

        # Fluxes and metabolite rates of change
        _v = pm.Deterministic("v", np.exp(_logkcat) * _f * _E, dims="reactions")
        _dm_dt = pm.Deterministic("dm_dt", _S @ _v, dims="metabolites")

        # Create observed variables for observed kcats, fixed fluxes, all metabolite dm_dts
        _logkcat_obs = pm.Normal("log_kcat_obs", mu=_logkcat[_observed_kcat_indices], sigma=0.1, observed=np.log(_observed_kcats))
        _v_obs = pm.Normal("v_obs", mu=_v[_observed_flux_indices], sigma=0.1, observed=_observed_fluxes)
        _dm_dt_obs = pm.Normal("dm_dt_obs", mu=_dm_dt, sigma=0.1, observed=np.zeros(4))

        # Sample
        idata_1 = pm.sample(draws=1000, chains=2)
    return (idata_1,)


@app.cell(hide_code=True)
def _(az, idata_1, plt):
    _axs = az.plot_trace(idata_1, combined=True)
    _fig = plt.gcf()

    _fig.set_size_inches(8, 6)
    _fig.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(f_prior_1, idata_1, np, plt):
    _fig, _ax = plt.subplots()

    _f_samples = idata_1.posterior["f"].data.flatten()

    _ax.hist(_f_samples, bins=50, density=True, label="Posterior")
    _ax.set_xlabel("$f$")
    _ax.set_ylabel("Density")

    _x = np.linspace(_f_samples.min(), _f_samples.max(), 100)
    _y = f_prior_1.pdf(_x)
    _ax_prior = _ax.twinx()
    _ax_prior.plot(_x, _y, color="tab:orange", label="Prior")
    _ax_prior.set_yticks([])

    _fig.legend(loc = "outside center right")
    _fig.suptitle("Prior-Posterior Plot for $f$")
    _fig.set_size_inches(4, 3)
    _fig.tight_layout()
    _fig
    return


@app.cell(hide_code=True)
def _(idata_1, kcat_marginal, logkcat_marginal, np, plt):
    # TODO: make these globals?
    _observed_kcats = np.array([10, 20])
    _observed_kcat_indices = np.array([3, 4])

    # Plot kcat prior-posterior plots
    _fig, (_logkcat_axs, _kcat_axs) = plt.subplots(2, 6)

    for _i, (_ax, _samples) in enumerate(
        zip(
            _logkcat_axs,
            idata_1.posterior["log_kcat"]
            .transpose("reactions", "chain", "draw")
            .data,
        )
    ):
        # Plot posterior histogram
        _samples = _samples.flatten()
        _ax.hist(_samples, bins=50, density=True)

        # Plot prior distribution
        _ax_prior = _ax.twinx()
        _x = np.linspace(_samples.min(), _samples.max(), 100)
        _y = logkcat_marginal.pdf(_x)
        _ax_prior.plot(_x, _y, color="tab:orange")

        # If there was an observed value, plot that
        if _i in _observed_kcat_indices:
            _ax_prior.vlines(np.log(_observed_kcats[_observed_kcat_indices == _i]), 0, _y.max(), color="r", linestyle="--")

        _ax.set_xlabel(rf"$\log k_{{cat}}^{{({_i})}}$")
        _ax_prior.set_yticks([])

    for _i, (_ax, _samples) in enumerate(
        zip(
            _kcat_axs,
            idata_1.posterior["kcat"]
            .transpose("reactions", "chain", "draw")
            .data,
        )
    ):
        # Plot posterior histogram
        _samples = _samples.flatten()
        _ax.hist(_samples, bins=50, density=True)

        # Plot prior distribution
        _ax_prior = _ax.twinx()
        _x = np.linspace(_samples.min(), _samples.max(), 100)
        _ax_prior.plot(_x, kcat_marginal.pdf(_x), color="tab:orange")

        # If there was an observed value, plot that
        if _i in _observed_kcat_indices:
            _ax_prior.vlines(_observed_kcats[_observed_kcat_indices == _i], 0, _y.max(), color="r", linestyle="--")

        _ax.set_xlabel(rf"$k_{{cat}}^{{({_i})}}$")
        _ax_prior.set_yticks([])

    _logkcat_axs[0].set_ylabel("Density")
    _kcat_axs[0].set_ylabel("Density")

    _fig.suptitle(r"Prior-Posterior Plots for $\log k_{cat}$, $k_{cat}$")
    _fig.set_size_inches(12, 5)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(az, idata_1):
    az.plot_pair(idata_1, var_names=["kcat", "f"], kind="scatter", marginals=True)
    return


if __name__ == "__main__":
    app.run()
