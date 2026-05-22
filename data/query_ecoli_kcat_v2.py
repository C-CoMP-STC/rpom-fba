#!/usr/bin/env python3
"""
query_ecoli_kcat.py  (v2)
=========================
Query kcat (turnover number) values measured in wild-type Escherichia coli
from two databases:

  • SABIO-RK  – REST API  (no credentials required)
  • BRENDA    – SOAP API  (free account required; register at brenda-enzymes.org)

Results are merged, metadata-enriched, and mapped to BiGG reaction IDs from
the iJO1366 E. coli K-12 genome-scale metabolic model for downstream FBA use.

Dependencies
------------
    pip install cobra          # COBRApy

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

    Additional CLI flags:
    --match-precedence  Comma-separated strategy order, default: uniprot,metabolites,ec
    --skip-participants  Skip SABIO-RK reaction participant fetch (disables metabolite matching)

Output columns (TSV)
--------------------
    source              SABIO-RK | BRENDA
    ec_number           EC number of the enzyme
    kcat_value          Measured turnover number (s⁻¹)
    kcat_end            Upper bound, if provided
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
    bigg_match_basis   Strategy used for BiGG mapping: uniprot | metabolites | ec | none
    substrate_mnx      MetaNetX IDs of SABIO-RK reaction substrates (semicolon-sep)
    product_mnx        MetaNetX IDs of SABIO-RK reaction products (semicolon-sep)
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import logging
import re
import time
from io import StringIO
from itertools import product as iproduct
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

import cobra
import pandas as pd
import requests

log = logging.getLogger(__name__)


# ============================================================
# SABIO-RK
# ============================================================

_SABIORK_ENTRY_IDS_URL = (
    "https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs"
)
_SABIORK_EXPORT_URL = "https://sabiork.h-its.org/entry/exportToExcelCustomizable"
_SABIORK_PARTICIPANTS_URL = (
    "https://sabiork.h-its.org/sabioRestWebServices/searchReactionParticipants"
)

# All fields we want back from SABIO-RK.
# "Parameter" expands to type/startValue/endValue/unit/associatedSpecies columns.
_SABIORK_FIELDS = [
    "EntryID",
    "Organism",
    "UniprotID",
    "ECNumber",
    "Parameter",  # kinetic parameter block; see _parse_sabiork_parameter()
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
        "Organism": f'"{organism}"',
        "parameterType": "kcat",
    }

    for attempt, extra in enumerate([{"EnzymeType": "wildtype"}, {}]):
        q_dict = {**base_query, **extra}
        q_str = " AND ".join(f"{k}:{v}" for k, v in q_dict.items())
        log.info("SABIO-RK entryID query: %s", q_str)

        resp = requests.get(
            _SABIORK_ENTRY_IDS_URL,
            params={"format": "txt", "q": q_str},
            timeout=60,
        )
        resp.raise_for_status()
        text = resp.text.strip()

        ids = [int(x) for x in text.splitlines() if x.strip().lstrip("-").isdigit()]
        if ids:
            log.info("  → %d entry IDs", len(ids))
            return ids

        if attempt == 0:
            log.warning(
                "EnzymeType:wildtype filter returned 0 entries; "
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
            i // batch_size + 1,
            n_batches,
            len(batch),
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
        time.sleep(0.4)  # be polite to the SABIO-RK servers

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
            "parameter.type": "param_type",
            "parameter.startValue": "kcat_value",
            "parameter.endValue": "kcat_value_end",
            "parameter.unit": "kcat_unit",
            "parameter.associatedSpecies": "substrate_param",
        }
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    # Layout B: single 'Parameter' column (pipe-delimited)
    elif "Parameter" in df.columns:
        parts = df["Parameter"].str.split("|", expand=True)
        col_names = [
            "param_type",
            "kcat_value",
            "kcat_value_end",
            "kcat_unit",
            "substrate_param",
        ]
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
        "ECNumber": "ec_number",
        "Organism": "organism",
        "UniprotID": "uniprot_id",
        "Substrate": "substrate",
        "pH": "pH",
        "Temperature": "temperature_C",
        "Buffer": "buffer",
        "CellularLocation": "cellular_location",
        "SabioReactionID": "sabiork_reaction_id",
        "KeggReactionID": "kegg_reaction_id",
        "ReactomeReactionID": "reactome_reaction_id",
        "PubMedID": "pubmed_id",
        "EntryID": "sabiork_entry_id",
        "kcat_value": "kcat_value",
        "kcat_value_end": "kcat_value_end",
        "kcat_unit": "kcat_unit",
        "substrate_param": "substrate_from_param",
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


def fetch_sabiork(
    organism: str = "Escherichia coli",
    cache_dir: Path = Path("."),
    use_cached: bool = False,
) -> pd.DataFrame:
    """
    Full SABIO-RK pipeline → canonical DataFrame.

    Results are cached to ``<cache_dir>/sabiork_df.tsv`` so
    subsequent runs can skip fetching (if use_cached option selected.)
    """
    cache_path = cache_dir / "sabiork_df.tsv"
    if use_cached:
        if cache_path.exists():
            df = pd.read_csv(cache_path, sep="\t")
            log.info("Read cached SABIO-RK kcat values.")
            return df
        else:
            log.warning("No cached kcat values from SABIO-RK found! Fetching...")
    
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

    # Save to cache
    df.to_csv(cache_path, sep="\t", index=False)
    log.info(f"Saved SABIO-RK results to {cache_path}.")

    return df


# ============================================================
# SABIO-RK reaction participants  (NEW in v2)
# ============================================================


def fetch_sabiork_reaction_participants(
    reaction_ids: List[int],
    cache_dir: Path = Path("."),
    sleep: float = 0.25,
) -> Dict[int, Dict[str, List[Dict[str, str]]]]:
    """
    For each unique SabioReactionID, retrieve the substrate and product
    compound names plus ChEBI / KEGG IDs via the SABIO-RK REST API.

    Results are cached to ``<cache_dir>/sabiork_participants_cache.json`` so
    subsequent runs skip already-fetched reactions.

    Returns
    -------
    dict mapping reaction_id (int) →
        {"substrates": [{"name": ..., "chebi": ..., "kegg": ...}, ...],
         "products":   [...]}
    """
    cache_path = cache_dir / "sabiork_participants_cache.json"
    cache: Dict[str, Any] = {}
    if cache_path.exists():
        with cache_path.open() as fh:
            cache = json.load(fh)

    to_fetch = [rid for rid in reaction_ids if str(rid) not in cache]
    log.info(
        "Fetching SABIO-RK participants: %d new reactions (%d already cached)",
        len(to_fetch),
        len(reaction_ids) - len(to_fetch),
    )

    for i, rid in enumerate(to_fetch):
        if i > 0 and i % 100 == 0:
            log.info("  … %d/%d participants fetched", i, len(to_fetch))
            with cache_path.open("w") as fh:
                json.dump(cache, fh)

        resp = requests.get(
            _SABIORK_PARTICIPANTS_URL,
            params={
                "SabioReactionID": str(rid),
                "fields[]": [
                    "Name",
                    "Role",
                    "ChebiID",
                    "KeggCompoundID",
                    "SabioCompoundID",
                ],
            },
            timeout=30,
        )

        entry: Dict[str, List] = {"substrates": [], "products": []}
        if resp.status_code == 200 and resp.text.strip():
            try:
                pdf = pd.read_csv(StringIO(resp.text), sep="\t", dtype=str)
                for _, row in pdf.iterrows():
                    role = str(row.get("Role", "")).strip().lower()
                    c = {
                        "name": str(row.get("Name", "")).strip(),
                        "chebi": str(row.get("ChebiID", "")).strip(),
                        "kegg": str(row.get("KeggCompoundID", "")).strip(),
                    }
                    # Roles include 'substrate', 'product',
                    # 'cofactor_substrate', 'cofactor_product', etc.
                    if "substrate" in role:
                        entry["substrates"].append(c)
                    elif "product" in role:
                        entry["products"].append(c)
            except Exception as exc:
                log.debug("Could not parse participants for reaction %d: %s", rid, exc)

        cache[str(rid)] = entry
        time.sleep(sleep)

    # Final cache write
    with cache_path.open("w") as fh:
        json.dump(cache, fh)

    return {rid: cache[str(rid)] for rid in reaction_ids if str(rid) in cache}


# ============================================================
# MetaNetX compound cross-reference  (NEW in v2)
# ============================================================

# MetaNetX chem_xref.tsv: columns are XREF, MNX_ID, DESCRIPTION (MNXref ≥ 4.0)
# XREF format: "chebi:1234"  (bare number, no 'CHEBI:' prefix)
#              "kegg.compound:C00002"
#              "bigg.metabolite:atp"
# Download page: https://www.metanetx.org/mnxdoc/mnxref.html
_METANETX_CHEM_XREF_URL = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv"

# IDs in iJO1366 may be outdated, so need to be mapped to latest MetaNetX IDs
# METANETX chem_depr.tsv: columns are deprecated_ID, ID, version
# (where `version` refers to the MNXref version in which the deprecated_ID was deprecated)
_METANETX_CHEM_DEPR_URL = "https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_depr.tsv"


def build_metanetx_chebi_map(
    cache_dir: Path = Path("."),
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Download (once) and parse MetaNetX ``chem_xref.tsv``.

    Returns
    -------
    chebi_to_mnx : dict  ``"CHEBI:1234"`` → ``"MNXM123"``
                   (keys normalised to upper-case CHEBI: prefixed form)
    kegg_to_mnx  : dict  ``"C00002"`` → ``"MNXM123"``
    """
    xref_path = cache_dir / "MetaNetX_e_chem_xref.tsv"
    if not xref_path.exists():
        log.info("Downloading MetaNetX chem_xref.tsv (~25 MB) …")
        resp = requests.get(_METANETX_CHEM_XREF_URL, timeout=300, stream=True)
        resp.raise_for_status()
        with xref_path.open("wb") as fh:
            for chunk in resp.iter_content(1 << 16):
                fh.write(chunk)
        log.info("  Cached → %s", xref_path)
    else:
        log.info("Using cached MetaNetX xref at %s", xref_path)

    chebi_to_mnx: Dict[str, str] = {}
    kegg_to_mnx: Dict[str, str] = {}

    with xref_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            xref, mnx_id = parts[0].strip(), parts[1].strip()
            if not mnx_id.startswith("MNX"):
                continue

            if xref.startswith("chebi:"):
                # MNXref format: "chebi:1234" (bare number).
                # We store as "CHEBI:1234" to match SABIO-RK's output.
                num = xref[6:]  # strip "chebi:"
                chebi_to_mnx[f"CHEBI:{num}"] = mnx_id
                chebi_to_mnx[num] = mnx_id  # also bare number fallback
            elif xref.startswith("kegg.compound:") or xref.startswith("kegg:"):
                kegg_id = xref.rsplit(":", 1)[-1]
                kegg_to_mnx[kegg_id] = mnx_id

    log.info(
        "MetaNetX: %d ChEBI → MNXM, %d KEGG → MNXM",
        len(chebi_to_mnx) // 2,  # divided by 2 because we store with/without prefix
        len(kegg_to_mnx),
    )
    return chebi_to_mnx, kegg_to_mnx


def enrich_sabiork_with_mnx(
    df: pd.DataFrame,
    participants: Dict[int, Dict],
    chebi_to_mnx: Dict[str, str],
    kegg_to_mnx: Dict[str, str],
) -> pd.DataFrame:
    """
    Add ``substrate_mnx`` and ``product_mnx`` columns to a SABIO-RK DataFrame.

    Each cell contains a semicolon-separated string of MetaNetX MNXM IDs for
    the measured reaction's substrates / products, resolved via the
    MetaNetX ChEBI/KEGG cross-reference.  Missing mappings produce empty strings.
    """

    def _compounds_to_mnx(compounds: List[Dict]) -> Set[str]:
        mnx_ids: Set[str] = set()
        for c in compounds:
            chebi = c.get("chebi", "").strip()
            kegg = c.get("kegg", "").strip()

            # Normalise SABIO-RK ChEBI IDs: accept "CHEBI:1234", "1234"
            if chebi and chebi not in ("nan", ""):
                # Ensure "CHEBI:" prefix
                if not chebi.upper().startswith("CHEBI:"):
                    chebi_norm = f"CHEBI:{chebi}"
                else:
                    chebi_norm = chebi.upper()
                if chebi_norm in chebi_to_mnx:
                    mnx_ids.add(chebi_to_mnx[chebi_norm])
                elif chebi_norm.split(":")[-1] in chebi_to_mnx:
                    mnx_ids.add(chebi_to_mnx[chebi_norm.split(":")[-1]])

            if kegg and kegg not in ("nan", "") and kegg in kegg_to_mnx:
                mnx_ids.add(kegg_to_mnx[kegg])

        return mnx_ids

    sub_mnx_col: List[str] = []
    prod_mnx_col: List[str] = []

    for _, row in df.iterrows():
        rid_raw = row.get("sabiork_reaction_id", None)
        sub_mnx: Set[str] = set()
        prod_mnx: Set[str] = set()
        if pd.notna(rid_raw):
            try:
                rid = int(float(rid_raw))
                p = participants.get(rid, {})
                sub_mnx = _compounds_to_mnx(p.get("substrates", []))
                prod_mnx = _compounds_to_mnx(p.get("products", []))
            except (ValueError, TypeError):
                pass
        sub_mnx_col.append(";".join(sorted(sub_mnx)))
        prod_mnx_col.append(";".join(sorted(prod_mnx)))

    df = df.copy()
    df["substrate_mnx"] = sub_mnx_col
    df["product_mnx"] = prod_mnx_col
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
    return (email, pw_hash)  # f"{email},{pw_hash}"


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
        "ec_number": fields.get("ecNumber"),
        "organism": fields.get("organism"),
        "kcat_value": fields.get("turnoverNumber"),
        "kcat_maximum": fields.get("turnoverNumberMaximum"),
        "kcat_unit": "1/s",
        "substrate": fields.get("substrate"),
        "pubmed_id": fields.get("literature"),
        "commentary": commentary,
        "textmining": fields.get("textmining"),
    }


def fetch_brenda(
    email: str, password: str, organism: str = "Escherichia coli"
) -> pd.DataFrame:
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
        log.warning(
            "BRENDA returned an empty result (check credentials / organism name)."
        )
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
# BiGG / iJO1366 mapping  (MODIFIED in v2)
# ============================================================

_IJO1366_URL = "http://bigg.ucsd.edu/static/models/iJO1366.xml"


# ---------------------------------------------------------------------------
# GPR catalyst extraction  (NEW helper, from user-supplied snippet)
# ---------------------------------------------------------------------------


def _catalysts_from_gpr(gpr: Any) -> Optional[List]:
    """
    Parse a COBRApy GPR object and return catalyst gene-ID lists.

    The return value is a list in which each element represents one
    iso-catalyst:
      - a plain ``str`` for a single-gene catalyst
      - a ``tuple`` of strs for a multi-gene complex  (AND relationship)

    Top-level OR nodes are handled; nested ORs raise ``NotImplementedError``
    (not present in iJO1366, per user confirmation).

    Returns ``None`` for reactions with no GPR (spontaneous / exchange).
    """

    def _catalysts_in(node: ast.expr) -> List:
        match type(node):
            case ast.BoolOp:
                if isinstance(node.op, ast.Or):
                    raise NotImplementedError(
                        "Nested OR in GPR — not expected in iJO1366"
                    )
                elif isinstance(node.op, ast.And):
                    return list(
                        iproduct(*(_catalysts_in(child) for child in node.values))
                    )
            case ast.Name:
                return [node.id]
            case _:
                raise ValueError(f"Unrecognised GPR node type: {type(node)}")

    root = gpr.body
    if root is None:
        return None

    if isinstance(root, ast.BoolOp) and isinstance(root.op, ast.Or):
        catalysts: List = []
        for child in root.values:
            catalysts += _catalysts_in(child)
        return catalysts
    else:
        return _catalysts_in(root)


# ---------------------------------------------------------------------------
# BiGG index dataclass  (NEW)
# ---------------------------------------------------------------------------


class BiggIndex:
    """
    Pre-built lookup structures derived from the iJO1366 model.

    Attributes
    ----------
    ec_to_rxns : dict  EC number → [rxn_id, ...]
    uniprot_to_rxns : dict  UniProt accession → [rxn_id, ...]
                      (any gene appearing in the reaction's GPR)
    sub_mnx_to_rxns : dict  MNXM_ID → {rxn_id, ...}   (inverted substrate index)
    prod_mnx_to_rxns : dict  MNXM_ID → {rxn_id, ...}   (inverted product index)
    rxn_details : dict  rxn_id → {
                      "ec_numbers":       set[str],
                      "sub_mnx":          frozenset[str],
                      "prod_mnx":         frozenset[str],
                      "reversible":       bool,
                      "catalyst_uniprots": list[frozenset[str]] | None,
                          # one frozenset per iso-catalyst; a single-gene
                          # catalyst is a frozenset of one UniProt ID;
                          # a complex is a frozenset of all subunit UniProt IDs
                  }
    """

    def __init__(self):
        self.ec_to_rxns: Dict[str, List[str]] = {}
        self.uniprot_to_rxns: Dict[str, List[str]] = {}
        self.sub_mnx_to_rxns: Dict[str, Set[str]] = {}
        self.prod_mnx_to_rxns: Dict[str, Set[str]] = {}
        self.rxn_details: Dict[str, Dict] = {}


def _update_model_mnx(model: cobra.Model, cache_dir: Path = Path(".")) -> cobra.Model:
    """
    Update the metabolite MetaNetX ids to the latest version.
    Downloads (once) and parses MetaNetX ``chem_depr.tsv``.
    """

    # Fetch chem_depr.tsv
    cache_path = cache_dir / "chem_depr.tsv"
    if not cache_path.exists():
        log.info("  Downloading MetaNetX chem_depr.tsv …")
        resp = requests.get(_METANETX_CHEM_DEPR_URL, timeout=300, stream=True)
        resp.raise_for_status()
        with cache_path.open("wb") as fh:
            for chunk in resp.iter_content(1 << 16):
                fh.write(chunk)
        log.info("    Cached → %s", cache_path)
    else:
        log.info("  Using cached MetaNetX chem_depr.tsv at %s", cache_path)

    # Parse into mapping
    mnx_old_to_new = {}
    with cache_path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            old_mnx, new_mnx = parts[0].strip(), parts[1].strip()
            mnx_old_to_new[old_mnx] = new_mnx

    # Update model
    n_updated_ids = 0
    for met in model.metabolites:
        mnx_old = met.annotation.get("metanetx.chemical")
        if mnx_old is None:
            continue
        
        mnx_new = mnx_old_to_new.get(mnx_old)
        if mnx_new is None:
            continue

        met.annotation["metanetx.chemical"] = mnx_new
        n_updated_ids += 1
    
    log.info("  Updated %d MetaNetX IDs to their latest versions.", n_updated_ids)
    return model


def _bigg_map_cobra(model_path: Path, cache_dir: Path = Path(".")) -> BiggIndex:
    """
    Build a ``BiggIndex`` from the iJO1366 SBML model using COBRApy.

    For each reaction:
    - EC numbers are read from ``rxn.annotation["ec-code"]``
    - Substrate/product MetaNetX IDs from ``met.annotation["metanetx.chemical"]``
    - Catalyst UniProt IDs by resolving GPR gene b-numbers via
      ``gene.annotation["uniprot"]`` (s0001 is skipped — empty annotation)
    """

    log.info("Parsing iJO1366 SBML with COBRApy …")
    model = cobra.io.read_sbml_model(str(model_path))

    # Update MetaNetX IDs to latest version
    model = _update_model_mnx(model, cache_dir)

    # Build gene-id → UniProt lookup first
    gene_id_to_uniprot: Dict[str, str] = {}
    for gene in model.genes:
        uniprot = gene.annotation.get("uniprot")
        if uniprot:
            # annotation value may be str or list
            uid = uniprot if isinstance(uniprot, str) else uniprot[0]
            gene_id_to_uniprot[gene.id] = uid

    index = BiggIndex()

    def _met_mnx(met) -> Set[str]:
        """Extract MetaNetX IDs from a metabolite annotation (str or list)."""
        raw = met.annotation.get("metanetx.chemical")
        if raw is None:
            return set()
        return {raw} if isinstance(raw, str) else set(raw)

    for rxn in model.reactions:
        rid = rxn.id

        # --- EC numbers ---
        ec_set: Set[str] = set()
        for ann_key in ("ec-code", "EC Number"):
            raw_ec = rxn.annotation.get(ann_key)
            if raw_ec is None:
                continue
            for ec in [raw_ec] if isinstance(raw_ec, str) else raw_ec:
                ec = ec.strip()
                if ec:
                    ec_set.add(ec)
                    index.ec_to_rxns.setdefault(ec, []).append(rid)

        # --- Metabolites ---
        sub_mnx: FrozenSet[str] = frozenset(
            mnx for met in rxn.reactants for mnx in _met_mnx(met)
        )
        prod_mnx: FrozenSet[str] = frozenset(
            mnx for met in rxn.products for mnx in _met_mnx(met)
        )
        for mnx in sub_mnx:
            index.sub_mnx_to_rxns.setdefault(mnx, set()).add(rid)
        for mnx in prod_mnx:
            index.prod_mnx_to_rxns.setdefault(mnx, set()).add(rid)

        # --- Catalyst UniProt sets ---
        catalyst_uniprots: Optional[List[FrozenSet[str]]] = None
        try:
            raw_cats = _catalysts_from_gpr(rxn.gpr)
            if raw_cats is not None:
                catalyst_uniprots = []
                for cat in raw_cats:
                    # cat is either a str (single gene) or a tuple (complex)
                    gene_ids = [cat] if isinstance(cat, str) else list(cat)
                    uniprots = frozenset(
                        gene_id_to_uniprot[g]
                        for g in gene_ids
                        if g in gene_id_to_uniprot
                    )
                    if uniprots:
                        catalyst_uniprots.append(uniprots)
        except (NotImplementedError, ValueError) as exc:
            log.debug("GPR parse issue for %s: %s", rid, exc)

        # --- Populate uniprot_to_rxns (any gene in GPR) ---
        for gene in rxn.genes:
            uid = gene_id_to_uniprot.get(gene.id)
            if uid:
                index.uniprot_to_rxns.setdefault(uid, []).append(rid)

        index.rxn_details[rid] = {
            "ec_numbers": ec_set,
            "sub_mnx": sub_mnx,
            "prod_mnx": prod_mnx,
            "reversible": rxn.lower_bound < 0,
            "catalyst_uniprots": catalyst_uniprots,
        }

    log.info(
        "  iJO1366 (COBRApy): %d reactions, %d EC entries, "
        "%d UniProt entries, %d substrate-MNXM entries",
        len(model.reactions),
        len(index.ec_to_rxns),
        len(index.uniprot_to_rxns),
        len(index.sub_mnx_to_rxns),
    )
    return index


def build_bigg_maps(cache_dir: Path = Path(".")) -> BiggIndex:
    """
    Download (once) and cache iJO1366, then build and return a ``BiggIndex``.
    [Return type changed in v2: was Tuple[Dict, Dict], now BiggIndex]
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

    return _bigg_map_cobra(model_path, cache_dir)


# ============================================================
# Multi-strategy BiGG ID assignment  (REWRITTEN in v2)
# ============================================================

_VALID_STRATEGIES = ("uniprot", "metabolites", "ec")
_DEFAULT_PRECEDENCE = ("ec", "uniprot", "metabolites")


def _match_uniprot(
    row: pd.Series, index: BiggIndex, candidates: Optional[List] = None
) -> List[str]:
    """Return BiGG reaction IDs whose GPR contains the row's UniProt ID."""
    uid = str(row.get("uniprot_id", "")).strip()
    if uid and uid != "nan":
        hits = index.uniprot_to_rxns.get(uid, [])
        return hits if candidates is None else sorted(set(candidates) & set(hits))
    return candidates if candidates is not None else []


def _match_metabolites(
    row: pd.Series, index: BiggIndex, candidates: Optional[List] = None
) -> List[str]:
    """
    Return BiGG reaction IDs whose substrate/product MetaNetX sets contain
    the row's measured substrates and products as subsets.

    Matching rule:
      sab_sub ⊆ rxn_sub  AND  sab_prod ⊆ rxn_prod
      (also checks the reversed direction for reversible reactions)
      At least one of sab_sub / sab_prod must be non-empty.

    Uses the inverted index (``sub_mnx_to_rxns``) for efficiency:
    the intersection of per-metabolite candidate sets is computed first,
    then the subset condition is verified.
    """

    def _parse_mnx(cell: Any) -> FrozenSet[str]:
        raw = str(cell) if pd.notna(cell) else ""
        ids = {x.strip() for x in raw.split(";") if x.strip() and x.strip() != "nan"}
        return frozenset(ids)

    sab_sub = _parse_mnx(row.get("substrate_mnx", ""))
    sab_prod = _parse_mnx(row.get("product_mnx", ""))

    if not sab_sub and not sab_prod:
        return []

    def _candidates_from_index(
        mnx_set: FrozenSet[str], inv: Dict[str, Set[str]]
    ) -> Optional[Set[str]]:
        """Intersect candidate sets for each MNXM in mnx_set."""
        cands: Optional[Set[str]] = None
        for mnx in mnx_set:
            rxns = inv.get(mnx, set())
            cands = rxns if cands is None else cands & rxns
        return cands  # None means mnx_set was empty

    def _check_direction(
        sub: FrozenSet, prod: FrozenSet, sub_inv: Dict, prod_inv: Dict
    ) -> Set[str]:
        cands_s = _candidates_from_index(sub, sub_inv)
        cands_p = _candidates_from_index(prod, prod_inv)

        # Intersect substrate and product candidate sets
        if cands_s is None and cands_p is None:
            return set()
        elif cands_s is None:
            merged = cands_p
        elif cands_p is None:
            merged = cands_s
        else:
            merged = cands_s & cands_p

        # Verify subset condition
        hits: Set[str] = set()
        for rid in merged or set():
            d = index.rxn_details.get(rid, {})
            sub_ok = (not sub) or sub.issubset(d.get("sub_mnx", frozenset()))
            prod_ok = (not prod) or prod.issubset(d.get("prod_mnx", frozenset()))
            if sub_ok and prod_ok:
                hits.add(rid)
        return hits

    # Forward direction
    hits = _check_direction(
        sab_sub,
        sab_prod,
        index.sub_mnx_to_rxns,
        index.prod_mnx_to_rxns,
    )

    # Reverse direction (substrates ↔ products) for reversible reactions
    if not hits and (sab_sub or sab_prod):
        rev_hits = _check_direction(
            sab_prod,
            sab_sub,
            index.sub_mnx_to_rxns,
            index.prod_mnx_to_rxns,
        )
        # Only accept reverse hits if the matched reaction is flagged reversible
        hits = {
            rid
            for rid in rev_hits
            if index.rxn_details.get(rid, {}).get("reversible", False)
        }

    return sorted(hits) if candidates is None else sorted(set(hits) & set(candidates))


def _match_ec(
    row: pd.Series, index: BiggIndex, candidates: Optional[List] = None
) -> List[str]:
    """Return BiGG reaction IDs sharing the row's EC number (exact + prefix)."""
    ec_raw = str(row.get("ec_number", "")).strip()
    if not ec_raw or ec_raw == "nan":
        return []
    if ec_raw in index.ec_to_rxns:
        hits = index.ec_to_rxns[ec_raw]
        return sorted(hits) if candidates is None else sorted(set(hits) & set(candidates))

    # Prefix match for partial ECs (e.g. "1.2.3.-")
    stem = ec_raw.rstrip("-").rstrip(".")
    hits: Set[str] = set()
    for stored_ec, rxn_ids in index.ec_to_rxns.items():
        if stored_ec.startswith(stem):
            hits.update(rxn_ids)
    return sorted(hits) if candidates is None else sorted(set(hits) & set(candidates))


_STRATEGY_FN = {
    "uniprot": _match_uniprot,
    "metabolites": _match_metabolites,
    "ec": _match_ec,
}


def add_bigg_ids(
    df: pd.DataFrame,
    index: BiggIndex,
    precedence: Tuple[str, ...] = _DEFAULT_PRECEDENCE,
) -> pd.DataFrame:
    """
    Add ``bigg_reaction_ids`` and ``bigg_match_basis`` columns.

    For each row the strategies in ``precedence`` are tried in order,
    filtering the results from the previous stage. The last strategy
    that returns > 0 matches is used, and the strategies tried are
    recorded in order in ``bigg_match_basis``.

    All matched BiGG reaction IDs are joined with semicolons.

    Parameters
    ----------
    precedence : tuple of strategy names in priority order.
        Valid values: "uniprot", "metabolites", "ec".
        Default: ("metabolites", "uniprot", "ec").
    """
    unknown = set(precedence) - set(_VALID_STRATEGIES)
    if unknown:
        raise ValueError(f"Unknown matching strategies: {unknown}")

    bigg_ids_col: List[str] = []
    match_basis_col: List[str] = []

    for _, row in df.iterrows():
        matched_ids: List[str] = []
        matched_basis: List[str] = []

        hits = None
        for strategy in precedence:
            fn = _STRATEGY_FN[strategy]
            next_hits = fn(row, index, candidates=hits)
            if len(next_hits) > 0:
                hits = next_hits
                matched_basis.append(strategy)
            else:
                break
        matched_ids = hits if hits is not None else []

        bigg_ids_col.append(";".join(sorted(set(matched_ids))))
        match_basis_col.append(
            ">".join(matched_basis) if len(matched_basis) > 0 else "none"
        )

    n_matched = sum(b != "none" for b in match_basis_col)
    basis_counts = pd.Series(match_basis_col).value_counts().to_dict()
    log.info(
        "BiGG mapping: %d/%d rows matched. Basis breakdown: %s",
        n_matched,
        len(df),
        basis_counts,
    )
    log.info("%d unique reactions mapped.", len(set(bigg_ids_col)))

    df = df.copy()
    df["bigg_reaction_ids"] = bigg_ids_col
    df["bigg_match_basis"] = match_basis_col
    return df


# ============================================================
# Output formatting  (MODIFIED: added bigg_match_basis, substrate_mnx, product_mnx)
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
    "bigg_match_basis",  # NEW
    "substrate_mnx",  # NEW
    "product_mnx",  # NEW
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
    for col in ("kcat_value", "kcat_value_end", "kcat_maximum", "pH", "temperature_C"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    present_lead = [c for c in _LEAD_COLS if c in df.columns]
    rest = [c for c in df.columns if c not in present_lead]
    df = df[present_lead + rest]

    df.to_csv(path, sep="\t", index=False)
    log.info("Saved %d rows → %s", len(df), path)


# ============================================================
# CLI  (MODIFIED: added --match-precedence, --skip-participants)
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
        "--brenda-email",
        default=None,
        help="BRENDA account e-mail. Required for BRENDA queries.",
    )
    p.add_argument(
        "--brenda-password",
        default=None,
        help="BRENDA account password.",
    )
    p.add_argument(
        "--organism",
        default="Escherichia coli",
        help="Organism name as it appears in BRENDA/SABIO-RK.",
    )
    p.add_argument("--use-cached", action="store_true", help="Use cached kcat values (if they exist).")
    p.add_argument("--skip-brenda", action="store_true", help="Skip BRENDA query.")
    p.add_argument("--skip-sabiork", action="store_true", help="Skip SABIO-RK query.")
    p.add_argument(
        "--skip-bigg",
        action="store_true",
        help="Skip BiGG ID mapping (avoids downloading iJO1366).",
    )
    p.add_argument(
        "--cache-dir",
        default="data/",
        help="Directory to cache stored data from previous runs.",
    )
    p.add_argument(
        "--output",
        default="data/ecoli_kcat.tsv",
        help="Output TSV file.",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    # NEW arguments
    p.add_argument(
        "--match-precedence",
        default=",".join(_DEFAULT_PRECEDENCE),
        help=(
            "Comma-separated matching strategy order. "
            "Strategies are tried left-to-right; the first with ≥1 hit wins. "
            f"Valid values: {', '.join(_VALID_STRATEGIES)}. "
            "Example: 'uniprot,metabolites,ec' (default) or 'metabolites,ec'."
        ),
    )
    p.add_argument(
        "--skip-participants",
        action="store_true",
        help=(
            "Skip the SABIO-RK reaction participant fetch. "
            "Disables 'metabolites' matching strategy (falls through to 'ec')."
        ),
    )
    return p


def main() -> None:
    args = _build_parser().parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%H:%M:%S",
    )

    # Parse and validate precedence
    precedence = tuple(s.strip() for s in args.match_precedence.split(",") if s.strip())
    unknown = set(precedence) - set(_VALID_STRATEGIES)
    if unknown:
        raise SystemExit(f"Unknown matching strategies: {unknown}")

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    frames: List[pd.DataFrame] = []

    # ── SABIO-RK ──────────────────────────────────────────
    sab: pd.DataFrame = pd.DataFrame()
    if not args.skip_sabiork:
        sab = fetch_sabiork(args.organism, cache_dir, args.use_cached)
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

    # ── Reaction participant fetch + MetaNetX enrichment ─
    # Only needed when metabolites strategy is active and we have SABIO-RK data
    if (
        not args.skip_bigg
        and not args.skip_participants
        and "metabolites" in precedence
        and not sab.empty
        and "sabiork_reaction_id" in sab.columns
    ):
        rid_series = (
            combined["sabiork_reaction_id"].dropna().apply(lambda x: int(float(x)))
        )
        unique_rids = sorted(rid_series.unique().tolist())
        log.info(
            "Fetching participants for %d unique SABIO-RK reactions …", len(unique_rids)
        )

        participants = fetch_sabiork_reaction_participants(unique_rids, cache_dir)
        chebi_to_mnx, kegg_to_mnx = build_metanetx_chebi_map(cache_dir)
        combined = enrich_sabiork_with_mnx(
            combined, participants, chebi_to_mnx, kegg_to_mnx
        )
    else:
        # Ensure columns exist even if empty (needed by _match_metabolites)
        if "substrate_mnx" not in combined.columns:
            combined["substrate_mnx"] = ""
        if "product_mnx" not in combined.columns:
            combined["product_mnx"] = ""

    # ── BiGG mapping ──────────────────────────────────────
    if not args.skip_bigg:
        # Since metabolite IDs in iJO1366 may be outdated, map to latest IDs

        bigg_index = build_bigg_maps(cache_dir)
        combined = add_bigg_ids(combined, bigg_index, precedence=precedence)
    else:
        combined["bigg_reaction_ids"] = ""
        combined["bigg_match_basis"] = "none"

    save_tsv(combined, args.output)
    log.info("Done.")


if __name__ == "__main__":
    main()
