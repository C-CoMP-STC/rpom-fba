"""Utility for accessing a local BioCyc R pom database.
An instance of Pathway Tools must be running with the -python option set, e.g. with

`> ./pathway-tools -lisp -python`
"""

import pythoncyc

ORGANISM = "|GCF_000011965|"


def access_rpom_biocyc(db=ORGANISM):
    try:
        return pythoncyc.select_organism(db)
    except Exception as e:
        raise ConnectionRefusedError(
            "Failed to connect to a Pathway Tools database. Did you run `pathway-tools -lisp -python?`") from e
