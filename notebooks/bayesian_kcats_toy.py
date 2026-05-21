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
    from scipy.stats import norm, lognorm, beta, gmean

    return az, beta, gmean, lognorm, norm, np, pd, plt, pm


@app.cell
def _(pd):
    # Load kcats, subset to E. coli wild-type
    kcats = pd.read_csv("data/ecoli_kcat.tsv", sep="\t")
    return (kcats,)


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
 
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
 
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Toy PyMC example

    Let's take the following simple bifurcating and merging reaction network:

    <img src="notebooks/public/toy_model.png" />

    We'll fix $v_1=10$, and observe $k_{cat}$'s for $v_4$ and $v_5$, leaving $k_{cat}$'s for $v_1, v_2, v_3,$ and $v_6$ unknown. We'll observe protein concentrations $E_i$ for all reactions. Ideally, as we vary $k_{cat}^{(4)}$ and $k_{cat}^{(5)}$, we should see $k_{cat}^{(2)}$ and $k_{cat}^{(3)}$ vary accordingly.

    Our stoichiometric matrix is:

    $$
    \begin{align*}
    S &= \begin{pmatrix}
        1 & -1 & -1 &    &    &   \\
          &  1 &    & -1 &    &   \\
          &    &  1 &    & -1 &   \\
          &    &    &  1 &  1 & -1
    \end{pmatrix}
    \end{align*}
    $$

    (where the omitted values $=0$).

    Each reaction will further have two associated random variables, $k_{cat}^{(i)}$ and $f_i$. For unobserved $k_{cat}$'s, we'll start with the uninformative prior $\log k_{cat}^{(i)} \sim \mathcal N(0.5, 2.5^2)$, informed by the marginal $k_{cat}$ distribution from above. For the $f_i$'s, we'll use $\text{Beta}$ priors. We'll first consider the case of a single shared $f$, then see how identifiability may be lost with individual $f_i$'s per reaction. Observed $k_{cat}$'s will be treated as deterministic.
    """)
    return


@app.cell
def _(gmean, kcats, np):
    kcats[~kcats["kcat_value"].isnull()].apply(
        lambda r: gmean(
            [r["kcat_value"]]
            + ([r["kcat_value_end"]] if not np.isnan(r["kcat_value_end"]) else [])
        ),
        axis=1,
    )
    return


@app.cell(hide_code=True)
def _(gmean, kcats, lognorm, norm, np, plt):
    # Fit marginal distribution for log-kcats
    _kcat_values = kcats[~kcats["kcat_value"].isnull()].apply(
        lambda r: gmean(
            [r["kcat_value"]]
            + ([r["kcat_value_end"]] if not np.isnan(r["kcat_value_end"]) else [])
        ),
        axis=1,
    )
    _kcat_values = _kcat_values[_kcat_values != 0]
    logkcat_marginal = norm(*norm.fit(np.log(_kcat_values)))
    kcat_marginal = lognorm(
        s=logkcat_marginal.std(), scale=np.exp(logkcat_marginal.mean())
    )

    _fig, (_ax1, _ax2) = plt.subplots(1, 2)

    _ax1.hist(np.log(_kcat_values), bins=50, density=True)
    _x = np.linspace(np.log(_kcat_values).min(), np.log(_kcat_values).max(), 100)
    _ax1.plot(_x, logkcat_marginal.pdf(_x))
    _ax1.set_xlabel(r"$\log k_{cat}$")
    _ax1.set_ylabel("Density")
    _ax1.set_title(
        rf"Marginal $\log k_{{cat}} \sim N ({logkcat_marginal.mean():.2f}, {logkcat_marginal.std():.2f})$"
    )

    _ax2.hist(_kcat_values, bins=50, density=True)
    _x = np.linspace(_kcat_values.min(), _kcat_values.max(), 100)
    _ax2.plot(_x, kcat_marginal.pdf(_x))
    _ax2.set_xlabel(r"$k_{cat}$")
    _ax2.set_yscale("log")
    _ax2.set_title(r"Marginal $k_{{cat}} \sim \log N(\cdot)$")

    _fig.set_size_inches(8, 2.5)
    _fig.tight_layout()
    _fig
    return kcat_marginal, logkcat_marginal


@app.cell(hide_code=True)
def _(mo):
    alpha_1 = mo.ui.slider(0.01, 10, 0.01, value=5, label=r"$\alpha$", show_value=True)
    beta_1 = mo.ui.slider(0.01, 10, 0.01, value=1, label=r"$\beta$", show_value=True)

    kcat_4_obs = mo.ui.number(1e-3, 1000, 1e-3, value=10, label="$k_{cat}^{(4),obs}$")
    kcat_5_obs = mo.ui.number(1e-3, 1000, 1e-3, value=20, label="$k_{cat}^{(5),obs}$")
    return alpha_1, beta_1, kcat_4_obs, kcat_5_obs


@app.cell(hide_code=True)
def _(alpha_1, beta, beta_1, kcat_4_obs, kcat_5_obs, mo, np, plt):
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

    mo.vstack([alpha_1, beta_1, _fig, kcat_4_obs, kcat_5_obs], heights="equal")
    return (f_prior_1,)


@app.cell(hide_code=True)
def _(alpha_1, beta_1, kcat_4_obs, kcat_5_obs, logkcat_marginal, np, pm):
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
    _observed_kcats = np.array([kcat_4_obs.value, kcat_5_obs.value])
    _observed_kcat_indices = np.array([3, 4])

    # Observed (fixed) fluxes
    _observed_fluxes = np.array([10])
    _observed_flux_indices = np.array([0])

    def build_model(S, E, observed_kcats, observed_kcat_idx, observed_fluxes, observed_flux_idx, sigma_log_kcat=0.1, sigma_v = 0.1, sigma_dm_dt = 0.1):
        _coords = {
            "reactions": np.arange(S.shape[1]),
            "metabolites": np.arange(S.shape[0]),
        }
        with pm.Model(coords=_coords) as model:
            # Priors for log-kcats
            _logkcat = pm.Normal("log_kcat", mu=logkcat_marginal.mean(), sigma=logkcat_marginal.std(), dims="reactions")
            _kcat = pm.Deterministic("kcat", np.exp(_logkcat), dims="reactions")

            # Prior for saturation f
            _f = pm.Beta("f", alpha=alpha_1.value, beta=beta_1.value)

            # Fluxes and metabolite rates of change
            _v = pm.Deterministic("v", np.exp(_logkcat) * _f * E, dims="reactions")
            _dm_dt = pm.Deterministic("dm_dt", S @ _v, dims="metabolites")

            # Create observed variables for observed kcats, fixed fluxes, all metabolite dm_dts
            _logkcat_obs = pm.Normal("log_kcat_obs", mu=_logkcat[observed_kcat_idx], sigma=sigma_log_kcat, observed=np.log(observed_kcats))
            _v_obs = pm.Normal("v_obs", mu=_v[observed_flux_idx], sigma=sigma_v, observed=observed_fluxes)
            _dm_dt_obs = pm.Normal("dm_dt_obs", mu=_dm_dt, sigma=sigma_dm_dt, observed=np.zeros(4))

        return model

    _model = build_model(_S, _E, _observed_kcats, _observed_kcat_indices, _observed_fluxes, _observed_flux_indices)
    idata_1 = pm.sample(model=_model, draws=1000, chains=2)
    return build_model, idata_1


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
def _(
    idata_1,
    kcat_4_obs,
    kcat_5_obs,
    kcat_marginal,
    logkcat_marginal,
    np,
    plt,
):
    _observed_kcats = np.array([kcat_4_obs.value, kcat_5_obs.value])
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


@app.cell(hide_code=True)
def _(build_model, kcat_4_obs, kcat_5_obs, kcat_marginal, np, plt, pm):
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
    _observed_kcats = np.array([kcat_4_obs.value, kcat_5_obs.value])
    _observed_kcat_indices = np.array([3, 4])

    # Observed (fixed) fluxes
    _observed_fluxes = np.array([10])
    _observed_flux_indices = np.array([0])



    # Plot kcat prior-posterior plots
    # _fig, _axs = plt.subplots(1, 7)
    # _kcat_axs = _axs[:6]
    # _lax = _axs[-1]

    _fig, _axs = plt.subplot_mosaic(
        [["k1", "k2", "k3", "lax"],
         ["k4", "k5", "k6", "lax"]]
    )
    _kcat_axs = [_axs[f"k{_i+1}"] for _i in range(6)]
    _lax = _axs["lax"]

    # Sample over varying noise in observed kcats
    _sigma_log_kcat = np.array([1e-2, 5e-2, 1e-1])
    _handles = []
    for _sigma in _sigma_log_kcat:

        # Build model with specified observation noise
        _model = build_model(_S, _E, _observed_kcats, _observed_kcat_indices, _observed_fluxes, _observed_flux_indices, sigma_log_kcat=_sigma)

        # Sample
        _idata = pm.sample(model=_model, draws=1000, chains=2)

        for _i, (_ax, _samples) in enumerate(
            zip(
                _kcat_axs,
                _idata.posterior["kcat"]
                .transpose("reactions", "chain", "draw")
                .data,
            )
        ):
            # Plot posterior histogram
            _samples = _samples.flatten()
            _handle = _ax.hist(_samples, bins=50, density=True, alpha=0.5, label=rf"$\sigma={_sigma:.1g}$")[2][0]

            # Build legend
            if _i == 0:
                _handles.append(_handle)


    # Plot priors
    for _i, (_ax, _samples) in enumerate(
            zip(
                _kcat_axs,
                _idata.posterior["kcat"]
                .transpose("reactions", "chain", "draw")
                .data,
            )
        ):
        # Plot prior distribution
        _ax_prior = _ax.twinx()
        _x = np.linspace(_samples.min(), _samples.max(), 100)
        _y = kcat_marginal.pdf(_x)
        _handle = _ax_prior.plot(_x, _y, color="tab:orange", label="Prior")[0]

        # If there was an observed value, plot that
        if _i in _observed_kcat_indices:
            _ax_prior.vlines(_observed_kcats[_observed_kcat_indices == _i], 0, _y.max(), color="r", linestyle="--")

        _ax.set_xlabel(rf"$k_{{cat}}^{{({_i+1})}}$")
        _ax_prior.set_yticks([])

        if _i == 0:
            _handles.append(_handle)

    _axs["k1"].set_ylabel("Density")
    _axs["k4"].set_ylabel("Density")

    # Plot legend
    _lax.legend(handles=_handles, loc="center left")
    _lax.axis("off")

    _fig.set_size_inches(8, 5)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(az, idata_1):
    az.plot_pair(idata_1, var_names=["kcat", "f"], kind="scatter", marginals=True)
    return


@app.cell(hide_code=True)
def _(idata_1, kcat_4_obs, kcat_5_obs, np, plt):
    import matplotlib as mpl
    from matplotlib.patches import Ellipse

    # Keep just the first chain
    _v = idata_1.posterior["v"][0, :, :].values.T

    # Coords of each metabolite
    _met_coords = {
        "_1": (-1, 0),
        "A": (0, 0),
        "B": (1, 1),
        "C": (1, -1),
        "X": (2, 0),
        "_2": (3, 0),
    }

    # Endpoints of each reaction
    _rxn_ends = {
        "v_1": ("_1", "A"),
        "v_2": ("A", "B"),
        "v_3": ("A", "C"),
        "v_4": ("B", "X"),
        "v_5": ("C", "X"),
        "v_6": ("X", "_2"),
    }

    _fig, (_ax1, _lax, _ax2) = plt.subplots(1, 3, gridspec_kw={"width_ratios": [1, 0.5, 2]})

    # Plot raw fluxes
    _handles = []
    for _i, _vi in enumerate(_v):
        _handles.append(_ax1.plot(_vi, label=f"$v_{_i+1}$")[0])
    _lax.legend(handles = _handles, loc="upper left")
    _lax.axis("off")

    _ax1.set_xlabel("Sample")
    _ax1.set_ylabel("Flux")

    # Set up colormap
    _cmap = mpl.colormaps["Reds"]
    _vmax = _v.max()
    _norm = mpl.colors.Normalize(0, _vmax)
    _normed_cmap = lambda x: _cmap(_norm(x))

    # Create nodes for metabolites
    for _met, _xy in _met_coords.items():
        if not _met.startswith("_"):
            _ax2.add_patch(Ellipse(_xy, 0.5, 0.5))
            _ax2.text(*_xy, _met, ha="center", va="center", color="w", size=20)

    # Draw reactions
    for (_rxn, (_a, _b)), _flux in zip(_rxn_ends.items(), _v):
        _from = _met_coords[_a]
        _to = _met_coords[_b]
        _midpoint = (np.array(_from) + np.array(_to)) / 2
        _text_offset_x = np.array([0.2, 0])
        _text_offset_y = np.array([0, 0.1])
        _left = _midpoint[0] <= 1
        _above_y0 = _midpoint[1] >= 0
        _ax2.annotate(
            "",
            _to,
            _from,
            arrowprops={
                "width": 5,
                "shrink": 0.2,
                "color": _normed_cmap(_flux.mean()),
            },
            zorder=-10,
        )
        _ax2.text(
            *(
                _midpoint
                + (-1 if _left else 1) * _text_offset_x
                + (1 if _above_y0 else -1) * _text_offset_y
            ),
            f"${_rxn}={_flux.mean():.2g}$",
            ha="center",
            va=("bottom" if _above_y0 else "top"),
        )

    _ax2.set_xlim(-1, 3)
    _ax2.set_ylim(-2, 2)
    _ax2.set_aspect("equal")
    _ax2.axis("off")

    _fig.suptitle(rf"$k_{{cat}}^{{(4), obs}} \leftarrow {kcat_4_obs.value:.2g}, k_{{cat}}^{{(5), obs}} \leftarrow {kcat_5_obs.value:.2g}$")

    _fig.colorbar(
        mpl.cm.ScalarMappable(norm=_norm, cmap=_cmap), ax=_ax2, label="Mean Flux"
    )

    _fig.set_size_inches(8, 4)
    _fig.tight_layout()
    _fig.subplots_adjust(wspace=0)
    _fig
    return


@app.cell
def _(build_model, np, pm):
    # Stoichiometric matrix
    _S = np.array([
        [1, -1, -1, 0, 0, 0],
        [0, 1, 0, -1, 0, 0],
        [0, 0, 1, 0, -1, 0],
        [0, 0, 0, 1, 1, -1]
    ])

    # Observed enzyme concentrations
    _E = np.array([1, 1, 1, 1, 1, 1])

    # Observed (fixed) fluxes
    _observed_fluxes = np.array([10])
    _observed_flux_indices = np.array([0])


    # Create grid
    _logkcat4 = np.linspace(-2, 2, 9)
    _logkcat5 = np.linspace(-2, 2, 9)
    _kcat4 = 10**(_logkcat4)
    _kcat5 = 10**(_logkcat5)

    F = np.zeros((_kcat4.size, _kcat5.size))

    for _i, _k4 in enumerate(_kcat4):
        for _j, _k5 in enumerate(_kcat5):
            _observed_kcats = np.array([_k4, _k5])
            _observed_kcat_indices = np.array([3, 4])

            _model = build_model(_S, _E, _observed_kcats, _observed_kcat_indices, _observed_fluxes, _observed_flux_indices)
            _idata = pm.sample(model=_model, draws=500, chains=1)

            F[_i, _j] = _idata.posterior["f"].mean()
    return (F,)


@app.cell
def _(F, np, plt):
    # Create grid
    _logkcat4 = np.linspace(-2, 2, 9)
    _logkcat5 = np.linspace(-2, 2, 9)
    _kcat4 = 10**_logkcat4
    _kcat5 = 10**_logkcat5

    _fig, _ax = plt.subplots()

    _im = _ax.imshow(
        F,
        extent=(
            _logkcat4.min(),
            _logkcat4.max(),
            _logkcat5.min(),
            _logkcat5.max(),
        ),
        origin="lower",
        vmin=0,
        vmax=1,
        cmap="Reds",
    )
    _ax.set_xticks(_logkcat4, [f"$10^{{{_k:.2g}}}$" for _k in _logkcat4])
    _ax.set_yticks(_logkcat4, [f"$10^{{{_k:.2g}}}$" for _k in _logkcat5])

    _ax.set_xlabel("$k_{cat}^{(4), obs}$")
    _ax.set_ylabel("$k_{cat}^{(5), obs}$")

    _fig.colorbar(_im, label="Saturation $f$")
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(build_model, np, plt, pm):
    # Plot v4, v5 as a function of k5 (k4 fixed at 10)

    # Stoichiometric matrix
    _S = np.array(
        [
            [1, -1, -1, 0, 0, 0],
            [0, 1, 0, -1, 0, 0],
            [0, 0, 1, 0, -1, 0],
            [0, 0, 0, 1, 1, -1],
        ]
    )

    # Observed enzyme concentrations
    _E = np.array([1, 1, 1, 1, 1, 1])

    # Observed (fixed) fluxes
    _observed_fluxes = np.array([10])
    _observed_flux_indices = np.array([0])

    _fig, _ax = plt.subplots()

    _logkcat5 = np.linspace(-2, 2, 9)
    _kcat5 = 10**_logkcat5
    _v4_map = []
    _v5_map = []
    for _x, _k5 in enumerate(_kcat5):
        _observed_kcats = np.array([10, _k5])
        _observed_kcat_indices = np.array([3, 4])

        _model = build_model(
            _S,
            _E,
            _observed_kcats,
            _observed_kcat_indices,
            _observed_fluxes,
            _observed_flux_indices,
        )
        _idata = pm.sample(model=_model, draws=500, chains=1)
        _map = pm.find_MAP(model=_model)

        # Store MAP
        _v4_map.append(_map["v"][3])
        _v5_map.append(_map["v"][4])

        # Plot samples
        _v4_samples = _idata.posterior["v"][:, :, 3].values.flatten()
        _ax.scatter(
            _x + np.random.normal(scale=0.01, size=_v4_samples.size),
            _v4_samples,
            color="0.8",
            alpha=0.01,
        )
        _v5_samples = _idata.posterior["v"][:, :, 4].values.flatten()
        _ax.scatter(
            _x + np.random.normal(scale=0.01, size=_v5_samples.size),
            _v5_samples,
            color="tab:blue",
            alpha=0.01,
        )

    _ax.plot(_v4_map, color="0.8", label="$v_4$")
    _ax.plot(_v5_map, label="$v_5$")

    _ax.legend()

    _ax.set_xticks(
        range(_logkcat5.size), [f"$10^{{{_k:.2g}}}$" for _k in _logkcat5]
    )
    _ax.set_xlabel("$k_{cat}^{(5), obs}$")
    _ax.set_ylabel("Flux")

    _fig.set_size_inches(6, 4)
    _fig.tight_layout()
    _fig
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
