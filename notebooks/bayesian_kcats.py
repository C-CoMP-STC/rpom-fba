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

    return load_model, pd


@app.cell
def _(pd):
    # Load kcats, subset to E. coli wild-type
    kcats = pd.read_csv("data/ecoli_kcat.tsv", sep="\t")

    # Load mapping of ecocyc gene ids to accessions, reactions
    ecocyc_genes_to_accessions = pd.read_csv("notebooks/data/EcoCyc_genes_to_reactions.tsv", sep="\t")
    ecocyc_genes_to_accessions["Reactions of gene"] = ecocyc_genes_to_accessions["Reactions of gene"].str.split(" // ")

    # Load Schmidt proteome dataset, merge in accessions
    prot = pd.read_csv("notebooks/data/schmidt2015_javier_table.tsv", sep="\t")
    prot = prot.merge(ecocyc_genes_to_accessions, left_on="EcoCycID", right_on="Gene Name")
    prot
    return (kcats,)


@app.cell
def _(kcats):
    kcats["bigg_reaction_ids"].value_counts()
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


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Plan:

    1. Check how well $f=1$ captures _E. coli_ dataset
        1. Take S
    """)
    return


@app.cell
def _(load_model):
    ecoli = load_model("iJO1366")

    return (ecoli,)


@app.cell
def _(ecoli, kcats):
    _ecoli_ids = [rxn.id for rxn in ecoli.reactions]
    _kcats_ids = kcats["bigg_reaction_ids"][~kcats["bigg_reaction_ids"].isnull()].unique()

    [_id for _id in _ecoli_ids if _id in _kcats_ids]
    return


@app.cell
def _(ecoli, kcats):
    _kcats_ids = kcats["bigg_reaction_ids"][~kcats["bigg_reaction_ids"].isnull()].unique()
    [rxn for rxn in ecoli.reactions if "and" in rxn.gene_reaction_rule and rxn.id in _kcats_ids]
    return


@app.cell
def _(ecoli):
    ecoli.genes.get_by_id('b2147').annotation
    return


@app.cell
def _(ecoli):
    [rxn for rxn in ecoli.reactions if len(rxn.genes) == 1][0]
    return


@app.cell
def _(ecoli):
    ecoli.reactions.get_by_id("14GLUCANtexi").gpr.body.id
    return


@app.cell
def _(product):
    list(product(["A", "B"], ["C"]))
    return


@app.cell
def _(ecoli):
    import ast
    from itertools import product

    def catalysts_from_gpr(gpr):
        def catalysts_in(node, catalysts=[]):
            match type(node):
                case ast.BoolOp:
                    if isinstance(node.op, ast.Or):
                        raise NotImplementedError()
                    elif isinstance(node.op, ast.And):
                        return list(product(*(catalysts_in(child) for child in node.values)))    
                case ast.Name:
                    return [node.id]
                case _:
                    raise ValueError(f"Unrecognized node type {type(node)}!")

        # Get root node
        root = gpr.body

        if root is None:
            return None

        # Top-level Or means several catalysts
        if isinstance(root, ast.BoolOp) and isinstance(root.op, ast.Or):
            catalysts = []
            for child in root.values:
                catalysts += catalysts_in(child)
            return catalysts
        else:
            return catalysts_in(gpr.body)

    catalysts = []
    for _rxn in ecoli.reactions:
        _catalysts = catalysts_from_gpr(_rxn.gpr)
        if _catalysts is not None:
            catalysts.append(_catalysts)
    catalysts
    return (product,)


@app.cell
def _(ecoli):
    from collections import Counter
    from functools import reduce
    from operator import add

    reduce(add, [Counter(_met.annotation.keys()) for _met in ecoli.metabolites])
    return Counter, add, reduce


@app.cell
def _(Counter, add, ecoli, reduce):
    len(ecoli.metabolites)

    reduce(add, [Counter(_gene.annotation.keys()) for _gene in ecoli.genes])
    return


@app.cell
def _(Counter, add, ecoli, reduce):
    reduce(add, [Counter(_met.annotation.keys()) for _met in ecoli.metabolites])
    return


@app.cell
def _(ecoli):
    ecoli.genes.get_by_id("s0001")
    return


@app.cell
def _(ecoli):
    [_gene for _gene in ecoli.genes if "uniprot" not in _gene.annotation]
    return


@app.cell
def _(ecoli):
    ecoli.reactions.get_by_id("CYSTRS")

    for _r in ecoli.metabolites.pyr_c.reactions:
        if ecoli.metabolites.aspsa_c in _r.metabolites:
            print(f"{_r}")

    # [m for m in ecoli.metabolites if "semialdehyde" in m.name.lower()]
    return


@app.cell
def _(ecoli):
    ecoli.reactions.DHDPS
    return


@app.cell
def _(ecoli):
    ecoli.genes.b2478.annotation
    return


@app.cell
def _(ecoli):
    [(met.name, met.annotation["metanetx.chemical"]) for met in ecoli.reactions.get_by_id("DHDPS").metabolites]
    return


@app.cell
def _(ecoli):
    ecoli.genes.get_by_id("b0526").annotation
    return


@app.cell
def _(Counter, ecoli):
    prefs = Counter()
    for _r in ecoli.reactions:
        _ec = _r.annotation.get("ec-code", [])
        if isinstance(_ec, str):
            _ec = [_ec]
        prefs += Counter(
            ".".join(_e.split(".")[:2])
            for _e in _ec
        )

    prefs
    return


if __name__ == "__main__":
    app.run()
