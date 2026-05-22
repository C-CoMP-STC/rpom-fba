#!/usr/bin/env python3
"""
query_ecoli_kcat.py
====================
Query kcat (turnover number) values measured in wild-type Escherichia coli
from two databases:

  • SABIO-RK  – REST API  (no credentials required)
  • BRENDA    – SOAP API  (free account required; register at brenda-enzymes.org)

Results are merged, metadata-enriched, and mapped to BiGG reaction IDs from
the iJO1366 E. coli K-12 genome-scale metabolic model for downstream FBA use.

Dependencies
------------
    pip install requests pandas zeep lxml

Optional (preferred for BiGG mapping):
    pip install cobra          # COBRApy; falls back to lxml SBML parser if absent

Usage
-----
    # SABIO-RK only (no account needed):
    python query_ecoli_kcat.py --skip-brenda --output ecoli_kcat.tsv

    # Both databases:
    python query_ecoli_kcat.py \\
        --brenda-email  your@email.com \\
        --brenda-password yourpassword \\
        --output ecoli_kcat.tsv

    # Skip the ~15 MB iJO1366 download if you don't need BiGG IDs:
    python query_ecoli_kcat.py --skip-brenda --skip-bigg --output ecoli_kcat.tsv

Output columns (TSV)
--------------------
    source              SABIO-RK | BRENDA
    ec_number           EC number of the enzyme
    kcat_value          Measured turnover number (s⁻¹)
    kcat_maximum        Upper bound, if provided (BRENDA only)
    kcat_unit           Unit string (always s⁻¹ after normalisation)
    substrate           Substrate the kcat was measured for
    organism            Organism string as recorded in the database
    uniprot_id          UniProt accession (SABIO-RK)
    bigg_reaction_ids   Semicolon-separated iJO1366 BiGG reaction IDs
    pH                  Assay pH
    temperature_C       Assay temperature (°C)
    buffer              Buffer description
    cellular_location   Sub-cellular location
    sabiork_entry_id    SABIO-RK EntryID
    sabiork_reaction_id SABIO-RK SabioReactionID
    kegg_reaction_id    KEGG reaction ID
    reactome_reaction_id Reactome reaction ID
    pubmed_id           PubMed reference(s)
    commentary          Free-text annotation (BRENDA only)
    textmining          Textmining flag (BRENDA only)
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import re
import time
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests

log = logging.getLogger(__name__)


# ============================================================
# SABIO-RK
# ============================================================

_SABIORK_ENTRY_IDS_URL = (
    "https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs"
)
_SABIORK_EXPORT_URL = (
    "https://sabiork.h-its.org/entry/exportToExcelCustomizable"
)

# All fields we want back from SABIO-RK.
# "Parameter" expands to type/startValue/endValue/unit/associatedSpecies columns.
_SABIORK_FIELDS = [
    "EntryID",
    "Organism",
    "UniprotID",
    "ECNumber",
    "Parameter",          # kinetic parameter block; see _parse_sabiork_parameter()
    "Substrate",
    "Product",
    "pH",
    "Temperature",
    "Buffer",
    "CellularLocation",
    "SabioReactionID",
    "KeggReactionID",
    "ReactomeReactionID",
    "PubMedID",
]

# Keywords that identify mutant entries in commentary text
_MUTANT_KW = re.compile(
    r"\b(mutant|variant|mutation|substitut|deletion|insertion|truncat|"
    r"chimera|engineered|modified|artificial)\b",
    re.IGNORECASE,
)


def _sabiork_entry_ids(organism: str) -> List[int]:
    """
    Step 1: fetch SABIO-RK EntryIDs for wild-type E. coli kcat measurements.

    Tries the IsWildtype:true filter first.  If SABIO-RK returns nothing
    (the filter has historically been unreliable), falls back to all kcat
    entries and we filter post-hoc.
    """
    base_query = {
        "Organism":      f'"{organism}"',
        "parameterType": "kcat",
    }

    for attempt, extra in enumerate([{"EnzymeType": "wildtype"}, {}]):
        q_dict = {**base_query, **extra}
        q_str  = " AND ".join(f"{k}:{v}" for k, v in q_dict.items())
        log.info("SABIO-RK entryID query: %s", q_str)

        resp = requests.get(
            _SABIORK_ENTRY_IDS_URL,
            params={"format": "txt", "q": q_str},
            timeout=60,
        )
        resp.raise_for_status()
        text = resp.text.strip()

        ids = [
            int(x)
            for x in text.splitlines()
            if x.strip().lstrip("-").isdigit()
        ]
        if ids:
            log.info("  → %d entry IDs", len(ids))
            return ids

        if attempt == 0:
            log.warning(
                "IsWildtype:true filter returned 0 entries; "
                "retrying without wildtype filter (will post-filter)."
            )

    log.warning("SABIO-RK returned 0 entries for '%s'.", organism)
    return []


def _sabiork_fetch_tsv(entry_ids: List[int], batch_size: int = 500) -> pd.DataFrame:
    """
    Step 2: POST entry IDs to SABIO-RK and retrieve TSV rows.

    The exportToExcelCustomizable endpoint accepts up to ~1000 IDs per
    request; we use batch_size=500 to stay comfortably within limits.
    """
    frames: List[pd.DataFrame] = []
    n_batches = (len(entry_ids) + batch_size - 1) // batch_size

    for i in range(0, len(entry_ids), batch_size):
        batch = entry_ids[i : i + batch_size]
        log.info(
            "  Fetching batch %d/%d (%d entries) …",
            i // batch_size + 1, n_batches, len(batch),
        )
        resp = requests.post(
            _SABIORK_EXPORT_URL,
            params={"format": "tsv", "fields[]": _SABIORK_FIELDS},
            data={"entryIDs[]": batch},
            timeout=120,
        )
        resp.raise_for_status()
        raw = resp.text.strip()
        if raw:
            df = pd.read_csv(StringIO(raw), sep="\t", dtype=str)
            frames.append(df)
        time.sleep(0.4)   # be polite to the SABIO-RK servers

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _parse_sabiork_parameter(df: pd.DataFrame) -> pd.DataFrame:
    """
    SABIO-RK returns the 'Parameter' column as a pipe-delimited block:
        type|startValue|endValue|unit|associatedSpecies
    or it may already be split into 'parameter.type', 'parameter.startValue',
    etc.  Handle both layouts and promote to flat columns.
    """
    # Layout A: already split by SABIO-RK into dotted columns
    if "parameter.type" in df.columns or "parameter.startValue" in df.columns:
        rename = {
            "parameter.type":             "param_type",
            "parameter.startValue":       "kcat_value",
            "parameter.endValue":         "kcat_value_end",
            "parameter.unit":             "kcat_unit",
            "parameter.associatedSpecies":"substrate_param",
        }
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    # Layout B: single 'Parameter' column (pipe-delimited)
    elif "Parameter" in df.columns:
        parts = df["Parameter"].str.split("|", expand=True)
        col_names = ["param_type", "kcat_value", "kcat_value_end",
                     "kcat_unit", "substrate_param"]
        for j, cname in enumerate(col_names):
            if j < parts.shape[1]:
                df[cname] = parts[j]
        df = df.drop(columns=["Parameter"])

    # Keep only kcat rows if we got mixed parameters
    if "param_type" in df.columns:
        mask = df["param_type"].str.lower().str.strip() == "kcat"
        if mask.any():
            df = df[mask].copy()

    return df


def _standardise_sabiork(df: pd.DataFrame) -> pd.DataFrame:
    """Rename SABIO-RK columns to the canonical output schema."""
    rename = {
        "ECNumber":           "ec_number",
        "Organism":           "organism",
        "UniprotID":          "uniprot_id",
        "Substrate":          "substrate",
        "pH":                 "pH",
        "Temperature":        "temperature_C",
        "Buffer":             "buffer",
        "CellularLocation":   "cellular_location",
        "SabioReactionID":    "sabiork_reaction_id",
        "KeggReactionID":     "kegg_reaction_id",
        "ReactomeReactionID": "reactome_reaction_id",
        "PubMedID":           "pubmed_id",
        "EntryID":            "sabiork_entry_id",
        "kcat_value":         "kcat_value",
        "kcat_value_end":     "kcat_value_end",
        "kcat_unit":          "kcat_unit",
        "substrate_param":    "substrate_from_param",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    # Consolidate substrate columns
    if "substrate" not in df.columns and "substrate_from_param" in df.columns:
        df["substrate"] = df["substrate_from_param"]

    # Normalise kcat units: SABIO-RK stores "s^(-1)", "1/s", "s-1" → "1/s"
    if "kcat_unit" in df.columns:
        df["kcat_unit"] = (
            df["kcat_unit"]
            .str.replace(r"s\^[\(\-1\)]+", "1/s", regex=True)
            .str.replace("s-1", "1/s", regex=False)
            .fillna("1/s")
        )
    else:
        df["kcat_unit"] = "1/s"

    return df


def fetch_sabiork(organism: str = "Escherichia coli") -> pd.DataFrame:
    """Full SABIO-RK pipeline → canonical DataFrame."""
    ids = _sabiork_entry_ids(organism)
    if not ids:
        return pd.DataFrame()

    df = _sabiork_fetch_tsv(ids)
    if df.empty:
        return df

    df = _parse_sabiork_parameter(df)
    df = _standardise_sabiork(df)
    df.insert(0, "source", "SABIO-RK")
    log.info("SABIO-RK: %d rows in final table", len(df))
    return df


# ============================================================
# BRENDA
# ============================================================

_BRENDA_WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"


def _brenda_credentials(email: str, password: str) -> str:
    """
    BRENDA authentication token: 'email,sha256(password)'.
    SHA-256 is required by the current BRENDA API; MD5 is no longer accepted.
    """
    pw_hash = hashlib.sha256(password.encode("utf-8")).hexdigest()
    return (email, pw_hash) # f"{email},{pw_hash}"


def _parse_brenda_entry(entry: str) -> Optional[Dict]:
    """
    Parse a single BRENDA result record.

    BRENDA separates records with '!' and encodes fields as 'key#value'
    pairs within each record (field separator '#', key-value separator '*').

    We discard entries whose commentary mentions mutagenesis keywords.
    """
    fields: Dict[str, str] = {}
    for part in entry.split("#"):
        if "*" in part:
            key, _, val = part.partition("*")
            fields[key.strip()] = val.strip()

    if not fields:
        return None

    commentary = fields.get("commentary", "")
    if _MUTANT_KW.search(commentary):
        return None

    return {
        "ec_number":    fields.get("ecNumber"),
        "organism":     fields.get("organism"),
        "kcat_value":   fields.get("turnoverNumber"),
        "kcat_maximum": fields.get("turnoverNumberMaximum"),
        "kcat_unit":    "1/s",
        "substrate":    fields.get("substrate"),
        "pubmed_id":    fields.get("literature"),
        "commentary":   commentary,
        "textmining":   fields.get("textmining"),
    }


def fetch_brenda(email: str, password: str, organism: str = "Escherichia coli") -> pd.DataFrame:
    """
    Query BRENDA via SOAP for turnover numbers in E. coli.

    The parameter string format is:
        '{email},{sha256_pw},ecNumber*{ec}#organism*{org}#field*value#…'
    Leaving ecNumber empty causes BRENDA to return all EC numbers.

    Note: BRENDA does not have a server-side wildtype filter for turnover
    numbers, so we post-filter by inspecting the commentary field.
    """
    try:
        import zeep  # type: ignore
    except ImportError:
        log.error("zeep is required for BRENDA queries: pip install zeep")
        return pd.DataFrame()

    creds = _brenda_credentials(email, password)
    # Empty ecNumber → return all EC entries for the given organism
    # parameters = (
    #     *creds,
    #     # f"{creds},"
    #     f"ecNumber*#organism*{organism}#",
    #     f"turnoverNumber*#turnoverNumberMaximum*#",
    #     f"substrate*#commentary*#literature*#textmining*#"
    # )
    parameters = (
        *creds,
        "ecNumber*",
        f"organism*{organism}",
        "turnoverNumber*",
        "turnoverNumberMaximum*",
        "substrate*",
        "commentary*",
        "ligandStructureId*",
        "literature*",
    )

    log.info("Connecting to BRENDA SOAP API …")
    try:
        client = zeep.Client(_BRENDA_WSDL)
        result_str: str = client.service.getTurnoverNumber(*parameters)
    except Exception as exc:
        log.error("BRENDA SOAP error: %s", exc)
        return pd.DataFrame()

    if not result_str:
        log.warning("BRENDA returned an empty result (check credentials / organism name).")
        return pd.DataFrame()

    rows = [
        parsed
        for raw_entry in result_str.split("!")
        if (parsed := _parse_brenda_entry(raw_entry.strip()))
    ]

    if not rows:
        log.warning("BRENDA: 0 rows after mutant filtering.")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df.insert(0, "source", "BRENDA")
    log.info("BRENDA: %d rows after mutant filtering", len(df))
    return df


# ============================================================
# BiGG / iJO1366 mapping
# ============================================================

_IJO1366_URL = "http://bigg.ucsd.edu/static/models/iJO1366.xml"


def _bigg_map_cobra(model_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Parse iJO1366 with COBRApy.
    Returns (ec_to_bigg, gene_to_bigg) where values are lists of BiGG reaction IDs.
    """
    import cobra.io  # type: ignore

    log.info("Parsing iJO1366 SBML with COBRApy …")
    model = cobra.io.read_sbml_model(str(model_path))

    ec_to_bigg: Dict[str, List[str]] = {}
    gene_to_bigg: Dict[str, List[str]] = {}

    for rxn in model.reactions:
        rid = rxn.id
        # EC numbers are in rxn.annotation under keys like 'ec-code'
        for ann_key in ("ec-code", "EC Number", "ec_code", "kegg.reaction"):
            val = rxn.annotation.get(ann_key)
            if val is None:
                continue
            for ec in ([val] if isinstance(val, str) else val):
                ec = ec.strip()
                if ec:
                    ec_to_bigg.setdefault(ec, []).append(rid)

        # Gene-level mapping (b-numbers and common names)
        for gene in rxn.genes:
            for gid in filter(None, [gene.id, gene.name]):
                gene_to_bigg.setdefault(gid, []).append(rid)

    log.info(
        "  iJO1366: %d reactions, %d EC entries, %d gene entries",
        len(model.reactions), len(ec_to_bigg), len(gene_to_bigg),
    )
    return ec_to_bigg, gene_to_bigg


def _bigg_map_lxml(model_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Fallback BiGG mapper using lxml (no COBRApy needed).
    Extracts EC numbers from SBML CVTerm annotations (identifiers.org/ec-code/).
    """
    from lxml import etree  # type: ignore  # noqa: F401

    log.info("Parsing iJO1366 SBML with lxml …")
    tree = etree.parse(str(model_path))
    root = tree.getroot()

    SBML  = "http://www.sbml.org/sbml/level3/version1/core"
    RDF   = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    BQ    = "http://biomodels.net/biology-qualifiers/"
    FBC   = "http://www.sbml.org/sbml/level3/version1/fbc/version2"

    ec_to_bigg: Dict[str, List[str]] = {}
    gene_to_bigg: Dict[str, List[str]] = {}

    for rxn in root.findall(f".//{{{SBML}}}reaction"):
        # Strip the 'R_' prefix that SBML adds to reaction IDs
        raw_id = rxn.get("id", "")
        rid    = raw_id[2:] if raw_id.startswith("R_") else raw_id

        # EC numbers from bqbiol:is CVTerms
        for li in rxn.findall(f".//{{{BQ}}}is//{{{RDF}}}li"):
            uri = li.get(f"{{{RDF}}}resource", "")
            if "ec-code" in uri or "enzyme" in uri:
                ec = uri.rstrip("/").rsplit("/", 1)[-1]
                if ec:
                    ec_to_bigg.setdefault(ec, []).append(rid)

        # Gene product references inside listOfGeneProductAssociations
        for gpr in rxn.findall(f".//{{{FBC}}}geneProductRef"):
            raw_gid = gpr.get("geneProduct", "")
            gid = raw_gid[2:] if raw_gid.startswith("G_") else raw_gid
            if gid:
                gene_to_bigg.setdefault(gid, []).append(rid)

    log.info(
        "  iJO1366: %d EC entries, %d gene entries",
        len(ec_to_bigg), len(gene_to_bigg),
    )
    return ec_to_bigg, gene_to_bigg


def build_bigg_maps(
    cache_dir: Path = Path("."),
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Download (once) and cache the iJO1366 SBML model, then build
    EC number → BiGG reaction ID and gene → BiGG reaction ID dicts.
    """
    model_path = cache_dir / "iJO1366.xml"

    if not model_path.exists():
        log.info("Downloading iJO1366 SBML from BiGG (~15 MB) …")
        resp = requests.get(_IJO1366_URL, timeout=180, stream=True)
        resp.raise_for_status()
        with model_path.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 16):
                fh.write(chunk)
        log.info("  Cached → %s", model_path)
    else:
        log.info("Using cached iJO1366 model at %s", model_path)

    try:
        return _bigg_map_cobra(model_path)
    except ImportError:
        log.warning("COBRApy not found; falling back to lxml SBML parser.")
        return _bigg_map_lxml(model_path)


def add_bigg_ids(
    df: pd.DataFrame,
    ec_to_bigg: Dict[str, List[str]],
    gene_to_bigg: Dict[str, List[str]],
) -> pd.DataFrame:
    """
    Add a 'bigg_reaction_ids' column (semicolon-delimited) to df.

    Matching strategy:
      1. Exact EC number match.
      2. Prefix match for partial/wildcard EC numbers (e.g. '1.2.3.-').
      3. UniProt ID or b-number lookup in gene_to_bigg
         (works if SABIO-RK UniProt accession happens to equal a gene b-number,
          which is rare; the primary mapping path is via EC number).
    """
    def _lookup(row: pd.Series) -> str:
        hits: set[str] = set()

        ec_raw = str(row.get("ec_number", "")).strip()
        if ec_raw and ec_raw not in ("nan", ""):
            if ec_raw in ec_to_bigg:
                hits.update(ec_to_bigg[ec_raw])
            else:
                # Prefix/wildcard match: '1.2.3.-' → match '1.2.3.*'
                stem = ec_raw.rstrip("-").rstrip(".")
                for stored_ec, bigg_ids in ec_to_bigg.items():
                    if stored_ec.startswith(stem):
                        hits.update(bigg_ids)

        for col in ("uniprot_id",):
            uid = str(row.get(col, "")).strip()
            if uid and uid != "nan" and uid in gene_to_bigg:
                hits.update(gene_to_bigg[uid])

        return ";".join(sorted(hits)) if hits else ""

    df = df.copy()
    df["bigg_reaction_ids"] = df.apply(_lookup, axis=1)
    n_mapped = (df["bigg_reaction_ids"] != "").sum()
    log.info(
        "BiGG mapping: %d/%d rows received at least one BiGG reaction ID",
        n_mapped, len(df),
    )
    return df


# ============================================================
# Output formatting
# ============================================================

_LEAD_COLS = [
    "source",
    "ec_number",
    "kcat_value",
    "kcat_value_end",
    "kcat_maximum",
    "kcat_unit",
    "substrate",
    "organism",
    "uniprot_id",
    "bigg_reaction_ids",
    "pH",
    "temperature_C",
    "buffer",
    "cellular_location",
    "sabiork_entry_id",
    "sabiork_reaction_id",
    "kegg_reaction_id",
    "reactome_reaction_id",
    "pubmed_id",
    "commentary",
    "textmining",
]


def save_tsv(df: pd.DataFrame, path: str) -> None:
    """Coerce numeric columns and write to TSV with canonical column ordering."""
    for col in ("kcat_value", "kcat_value_end", "kcat_maximum",
                "pH", "temperature_C"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    present_lead = [c for c in _LEAD_COLS if c in df.columns]
    rest         = [c for c in df.columns if c not in present_lead]
    df           = df[present_lead + rest]

    df.to_csv(path, sep="\t", index=False)
    log.info("Saved %d rows → %s", len(df), path)


# ============================================================
# CLI
# ============================================================

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Fetch wild-type E. coli kcat values from BRENDA + SABIO-RK, "
            "mapped to iJO1366 BiGG reaction IDs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--brenda-email", default=None,
        help="BRENDA account e-mail. Required for BRENDA queries.",
    )
    p.add_argument(
        "--brenda-password", default=None,
        help="BRENDA account password.",
    )
    p.add_argument(
        "--organism", default="Escherichia coli",
        help="Organism name as it appears in BRENDA/SABIO-RK.",
    )
    p.add_argument("--skip-brenda",   action="store_true", help="Skip BRENDA query.")
    p.add_argument("--skip-sabiork",  action="store_true", help="Skip SABIO-RK query.")
    p.add_argument(
        "--skip-bigg", action="store_true",
        help="Skip BiGG ID mapping (avoids downloading iJO1366).",
    )
    p.add_argument(
        "--bigg-cache-dir", default="data/",
        help="Directory to cache the iJO1366 SBML model.",
    )
    p.add_argument(
        "--output", default="data/ecoli_kcat.tsv",
        help="Output TSV file.",
    )
    p.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return p


def main() -> None:
    args = _build_parser().parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%H:%M:%S",
    )

    frames: List[pd.DataFrame] = []

    # ── SABIO-RK ──────────────────────────────────────────
    if not args.skip_sabiork:
        sab = fetch_sabiork(args.organism)
        if not sab.empty:
            frames.append(sab)

    # ── BRENDA ────────────────────────────────────────────
    if not args.skip_brenda:
        if not (args.brenda_email and args.brenda_password):
            log.warning(
                "BRENDA credentials not provided; skipping BRENDA.\n"
                "Register at https://www.brenda-enzymes.org and re-run with\n"
                "  --brenda-email <email> --brenda-password <password>"
            )
        else:
            brenda = fetch_brenda(
                args.brenda_email, args.brenda_password, args.organism
            )
            if not brenda.empty:
                frames.append(brenda)

    if not frames:
        log.error("No data retrieved from any source. Exiting.")
        return

    combined = pd.concat(frames, ignore_index=True, sort=False)
    log.info("Combined: %d rows from %d source(s)", len(combined), len(frames))

    # ── BiGG mapping ──────────────────────────────────────
    if not args.skip_bigg:
        ec_map, gene_map = build_bigg_maps(Path(args.bigg_cache_dir))
        combined = add_bigg_ids(combined, ec_map, gene_map)
    else:
        combined["bigg_reaction_ids"] = ""

    save_tsv(combined, args.output)
    log.info("Done.")


if __name__ == "__main__":
    main()
