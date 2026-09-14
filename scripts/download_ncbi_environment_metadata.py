#!/usr/bin/env python3

"""
Download NCBI BioSample environmental metadata for isolate genomes
in the mOTUs-db v4.0 thallusin analysis.

This version prioritizes speed:
- parallel requests with ThreadPoolExecutor
- no global request throttling
- retries only when NCBI returns 429 or temporary server errors
- checkpoint/resume support

Input
-----
output_tables/
    thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster.xlsx

Output
------
input_tables/
    mOTUs_GENO_NCBI_biosample_environment.tsv

Checkpoint
----------
.cache/
    ncbi_environment_download_status.tsv

Run from repository root:

    python scripts/download_ncbi_environment_metadata.py
"""

from __future__ import annotations

import os
import random
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_GENOME_TABLE = (
    PROJECT_ROOT
    / "output_tables"
    / "thallusin_mOTUsDB_124293genomes_mmseqs2_thal_cluster.xlsx"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "input_tables"
    / "mOTUs_GENO_NCBI_biosample_environment.tsv"
)

CACHE_DIR = PROJECT_ROOT / ".cache"

STATUS_FILE = (
    CACHE_DIR
    / "ncbi_environment_download_status.tsv"
)


# ============================================================
# Download settings
# ============================================================

NCBI_API_KEY = os.getenv("NCBI_API_KEY")

MAX_WORKERS = 3
REQUEST_TIMEOUT = 30
MAX_RETRIES = 8

CHECKPOINT_EVERY = 100


# ============================================================
# Missing-value handling
# ============================================================

MISSING_VALUES = {
    "",
    "missing",
    "not applicable",
    "not collected",
    "not provided",
    "not available",
    "not available: to be reported later",
    "unknown",
    "n/a",
    "na",
    "none",
    "not provided; submitted under migs 2.1",
}


# ============================================================
# Identifier conversion
# ============================================================

def genome_to_accession(genome: str) -> str | None:
    """
    Convert an mOTUs GENO identifier into an NCBI assembly accession.

    Example
    -------
    RSGB23-1_GCF-000023865-V1_GENO_10000001
    -> GCF_000023865.1
    """

    match = re.search(
        r"(GCF|GCA)-(\d+)-V(\d+)",
        genome,
    )

    if match is None:
        return None

    prefix, number, version = match.groups()

    return f"{prefix}_{number}.{version}"


# ============================================================
# Attribute normalization
# ============================================================

def normalize_attribute_name(name: str) -> str:
    """
    Normalize BioSample attribute names.
    """

    name = str(name).strip().lower()

    name = re.sub(
        r"[\s\-/]+",
        "_",
        name,
    )

    name = re.sub(
        r"[()]",
        "",
        name,
    )

    name = re.sub(
        r"_+",
        "_",
        name,
    )

    return name.strip("_")

# NCBI BioSample uses several synonymous or differently formatted
# field names for similar environmental metadata. These are mapped
# to a small canonical vocabulary so that downstream classification
# can use one consistent field name per type of evidence while
# retaining the original NCBI field name in `attribute_name`.

def canonical_attribute(name: str) -> str | None:
    """
    Map synonymous BioSample fields onto the compact
    environmental vocabulary used downstream.
    """

    n = normalize_attribute_name(name)

    mapping = {

        "isolation_source":
            "isolation_source",

        "direct_isolation_source":
            "isolation_source",

        "isolation_site":
            "isolation_site",

        "environment":
            "environment",

        "habitat":
            "habitat",

        "sample_habitat":
            "habitat",

        "env_broad_scale":
            "env_broad_scale",

        "env_broad_scale_biome":
            "env_broad_scale",

        "broad_scale_environmental_context":
            "env_broad_scale",

        "env_local_scale":
            "env_local_scale",

        "env_local_scale_feature":
            "env_local_scale",

        "local_environmental_context":
            "env_local_scale",

        "env_medium":
            "env_medium",

        "environmental_medium":
            "env_medium",

        "host":
            "host",

        "direct_host":
            "host",

        "host_body_habitat":
            "host_body_habitat",

        "body_habitat":
            "host_body_habitat",

        "plant_body_site":
            "plant_body_site",

        "sample_type":
            "sample_type",

        "source_type":
            "source_type",

        "env_package":
            "env_package",

        "isolation_location":
            "isolation_location",
    }

    if n in mapping:
        return mapping[n]

    environmental_packages = {
        "soil_environmental_package",
        "water_environmental_package",
        "host_associated_environmental_package",
        "wastewater_sludge_environmental_package",
        "sediment_environmental_package",
        "plant_associated_environmental_package",
    }

    if n in environmental_packages:
        return "environment_package"

    return None


def is_missing_value(value: str) -> bool:
    """
    Return True if a metadata value should be treated as missing.
    """

    return (
        str(value)
        .strip()
        .lower()
        in MISSING_VALUES
    )


# ============================================================
# HTTP helper
# ============================================================

def get_headers() -> dict:
    """
    Build request headers.

    If NCBI_API_KEY is set in the environment, it is used
    automatically without storing it in the repository.
    """

    headers = {
        "User-Agent":
            "thallusin-motus-analysis/1.0"
    }

    if NCBI_API_KEY:
        headers["api-key"] = NCBI_API_KEY

    return headers


# ============================================================
# Query NCBI
# ============================================================

def fetch_dataset_report(
    genome: str,
    accession: str,
) -> tuple[list[dict], str]:
    """
    Retrieve environmental BioSample metadata for one assembly.

    Returns
    -------
    rows
        Environmental metadata records.

    status
        "done" if the NCBI query succeeded, including cases where
        no useful environmental metadata were available.

        "failed" only after all retries were exhausted.
    """

    url = (
        "https://api.ncbi.nlm.nih.gov/"
        "datasets/v2alpha/genome/accession/"
        f"{accession}/dataset_report"
    )

    headers = get_headers()

    data = None


    for attempt in range(
        1,
        MAX_RETRIES + 1,
    ):

        try:

            response = requests.get(
                url,
                headers=headers,
                timeout=REQUEST_TIMEOUT,
            )


            # ------------------------------------------------
            # NCBI rate limit
            # ------------------------------------------------

            if response.status_code == 429:

                retry_after = response.headers.get(
                    "Retry-After"
                )

                if retry_after is not None:

                    try:
                        wait_seconds = float(
                            retry_after
                        )

                    except ValueError:
                        wait_seconds = 2 ** attempt

                else:

                    wait_seconds = min(
                        60,
                        2 ** attempt,
                    )


                wait_seconds += random.uniform(
                    0,
                    1,
                )


                print(
                    f"RATE LIMIT {accession}: "
                    f"waiting {wait_seconds:.1f} s "
                    f"(attempt {attempt}/{MAX_RETRIES})"
                )

                time.sleep(
                    wait_seconds
                )

                continue


            # ------------------------------------------------
            # Temporary server errors
            # ------------------------------------------------

            if response.status_code in {
                500,
                502,
                503,
                504,
            }:

                wait_seconds = min(
                    60,
                    2 ** attempt,
                )

                wait_seconds += random.uniform(
                    0,
                    1,
                )

                print(
                    f"SERVER ERROR "
                    f"{response.status_code} "
                    f"{accession}: "
                    f"waiting {wait_seconds:.1f} s"
                )

                time.sleep(
                    wait_seconds
                )

                continue


            response.raise_for_status()

            data = response.json()

            break


        except (
            requests.exceptions.Timeout,
            requests.exceptions.ConnectionError,
        ) as exc:

            wait_seconds = min(
                60,
                2 ** attempt,
            )

            wait_seconds += random.uniform(
                0,
                1,
            )

            print(
                f"NETWORK ERROR {accession}: "
                f"{exc}; waiting {wait_seconds:.1f} s"
            )

            time.sleep(
                wait_seconds
            )


        except Exception as exc:

            if attempt == MAX_RETRIES:

                print(
                    f"FAILED {accession}: {exc}"
                )

                return [], "failed"

            wait_seconds = min(
                60,
                2 ** attempt,
            )

            time.sleep(
                wait_seconds
            )


    if data is None:

        print(
            f"FAILED {accession}: "
            "maximum retries reached"
        )

        return [], "failed"


    # ========================================================
    # Parse report
    # ========================================================

    reports = data.get(
        "reports",
        [],
    )

    if not reports:
        return [], "done"


    report = reports[0]

    assembly_info = report.get(
        "assembly_info",
        {},
    )

    biosample = assembly_info.get(
        "biosample",
        {},
    )

    if not biosample:
        return [], "done"


    biosample_accession = biosample.get(
        "accession"
    )

    rows = []


    # ========================================================
    # General BioSample attributes
    # ========================================================

    for attribute in biosample.get(
        "attributes",
        [],
    ):

        name = attribute.get(
            "name"
        )

        value = attribute.get(
            "value"
        )

        if name is None or value is None:
            continue

        canonical = canonical_attribute(
            name
        )

        if canonical is None:
            continue

        value_clean = str(
            value
        ).strip()

        if is_missing_value(
            value_clean
        ):
            continue

        rows.append(
            {
                "genome":
                    genome,

                "assembly_accession":
                    accession,

                "biosample":
                    biosample_accession,

                "canonical_attribute":
                    canonical,

                "attribute_name":
                    name,

                "attribute_value":
                    value_clean,
            }
        )


    # ========================================================
    # Direct BioSample fields
    # ========================================================

    direct_fields = {

        "isolation_source":
            "isolation_source",

        "host":
            "host",
    }


    for field, canonical in direct_fields.items():

        value = biosample.get(
            field
        )

        if value is None:
            continue

        value_clean = str(
            value
        ).strip()

        if is_missing_value(
            value_clean
        ):
            continue

        rows.append(
            {
                "genome":
                    genome,

                "assembly_accession":
                    accession,

                "biosample":
                    biosample_accession,

                "canonical_attribute":
                    canonical,

                "attribute_name":
                    f"direct_{field}",

                "attribute_value":
                    value_clean,
            }
        )


    return rows, "done"


# ============================================================
# Save functions
# ============================================================

def save_results(rows: list[dict]) -> None:
    """
    Save metadata output atomically.
    """

    columns = [
        "genome",
        "assembly_accession",
        "biosample",
        "canonical_attribute",
        "attribute_name",
        "attribute_value",
    ]


    if rows:

        result = pd.DataFrame(
            rows
        )

        result = (
            result
            .drop_duplicates(
                subset=[
                    "genome",
                    "canonical_attribute",
                    "attribute_value",
                ]
            )
            .sort_values(
                [
                    "genome",
                    "canonical_attribute",
                    "attribute_value",
                ]
            )
            .reset_index(
                drop=True
            )
        )

    else:

        result = pd.DataFrame(
            columns=columns
        )


    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


    tmp_file = OUTPUT_FILE.with_suffix(
        ".tmp.tsv"
    )


    result.to_csv(
        tmp_file,
        sep="\t",
        index=False,
    )


    tmp_file.replace(
        OUTPUT_FILE
    )


def save_status(status_rows: list[dict]) -> None:
    """
    Save checkpoint file atomically.
    """

    CACHE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


    status = pd.DataFrame(
        status_rows
    )


    if not status.empty:

        status = (
            status
            .drop_duplicates(
                subset=["genome"],
                keep="last",
            )
            .sort_values(
                "genome"
            )
        )


    tmp_file = STATUS_FILE.with_suffix(
        ".tmp.tsv"
    )


    status.to_csv(
        tmp_file,
        sep="\t",
        index=False,
    )


    tmp_file.replace(
        STATUS_FILE
    )


# ============================================================
# Main
# ============================================================

def main() -> None:

    # --------------------------------------------------------
    # Check input
    # --------------------------------------------------------

    if not INPUT_GENOME_TABLE.exists():

        raise FileNotFoundError(
            "\nRequired input file not found:\n"
            f"{INPUT_GENOME_TABLE}\n\n"
            "Run the main analysis first:\n\n"
            "    Rscript scripts/thallusin_analysis.R\n"
        )


    print(
        f"Reading: {INPUT_GENOME_TABLE}"
    )

    print()

    if NCBI_API_KEY:
        print(
            "NCBI API key detected."
        )
    else:
        print(
            "No NCBI API key detected."
        )

    print(
        f"Parallel workers: {MAX_WORKERS}"
    )

    print()


    # --------------------------------------------------------
    # Read genome identifiers
    # --------------------------------------------------------

    genomes = pd.read_excel(
        INPUT_GENOME_TABLE,
        usecols=["genome"],
    )


    # --------------------------------------------------------
    # Select GENO genomes
    # --------------------------------------------------------

    genomes = (
        genomes.loc[
            genomes["genome"]
            .astype(str)
            .str.contains(
                "_GENO_",
                na=False,
            ),
            ["genome"],
        ]
        .drop_duplicates()
        .copy()
    )


    print(
        "GENO genomes:",
        f"{len(genomes):,}",
    )


    # --------------------------------------------------------
    # Parse assembly accessions
    # --------------------------------------------------------

    genomes[
        "assembly_accession"
    ] = genomes[
        "genome"
    ].map(
        genome_to_accession
    )


    missing_accessions = (
        genomes[
            "assembly_accession"
        ]
        .isna()
        .sum()
    )


    if missing_accessions:

        print(
            "GENO IDs without parsable NCBI accessions:",
            f"{missing_accessions:,}",
        )


    genomes = genomes.dropna(
        subset=[
            "assembly_accession"
        ]
    )


    print(
        "GENO genomes with NCBI accessions:",
        f"{len(genomes):,}",
    )


    # ========================================================
    # Resume previous run
    # ========================================================

    previous_rows = []
    status_rows = []


    if OUTPUT_FILE.exists():

        try:

            previous = pd.read_csv(
                OUTPUT_FILE,
                sep="\t",
            )

            previous_rows = previous.to_dict(
                orient="records"
            )

            print(
                "Existing metadata rows:",
                f"{len(previous_rows):,}",
            )

        except Exception as exc:

            print(
                "Could not read previous output:",
                exc,
            )


    completed_genomes = set()


    if STATUS_FILE.exists():

        try:

            previous_status = pd.read_csv(
                STATUS_FILE,
                sep="\t",
            )


            status_rows = previous_status.to_dict(
                orient="records"
            )


            completed_genomes = set(
                previous_status.loc[
                    previous_status["status"]
                    == "done",
                    "genome",
                ]
            )


            print(
                "Completed genomes from checkpoint:",
                f"{len(completed_genomes):,}",
            )

        except Exception as exc:

            print(
                "Could not read checkpoint:",
                exc,
            )


    # --------------------------------------------------------
    # Remaining genomes
    # --------------------------------------------------------

    remaining = genomes.loc[
        ~genomes["genome"].isin(
            completed_genomes
        )
    ].copy()


    print(
        "Assemblies remaining:",
        f"{len(remaining):,}",
    )


    if len(remaining) == 0:

        print(
            "\nNothing to download."
        )

        return


    # ========================================================
    # Download
    # ========================================================

    all_rows = list(
        previous_rows
    )


    with ThreadPoolExecutor(
        max_workers=MAX_WORKERS
    ) as executor:

        futures = {

            executor.submit(
                fetch_dataset_report,
                row.genome,
                row.assembly_accession,
            ):
                (
                    row.genome,
                    row.assembly_accession,
                )

            for row in remaining.itertuples(
                index=False
            )
        }


        total = len(
            futures
        )

        failed = 0


        for i, future in enumerate(
            as_completed(
                futures
            ),
            start=1,
        ):

            genome, accession = futures[
                future
            ]


            try:

                rows, status = future.result()

            except Exception as exc:

                print(
                    f"UNEXPECTED FAILURE "
                    f"{accession}: {exc}"
                )

                rows = []
                status = "failed"


            all_rows.extend(
                rows
            )


            status_rows.append(
                {
                    "genome":
                        genome,

                    "assembly_accession":
                        accession,

                    "status":
                        status,

                    "metadata_rows":
                        len(rows),
                }
            )


            if status == "failed":
                failed += 1


            # ------------------------------------------------
            # Checkpoint
            # ------------------------------------------------

            if (
                i % CHECKPOINT_EVERY == 0
                or i == total
            ):

                save_results(
                    all_rows
                )

                save_status(
                    status_rows
                )


                print(
                    f"Processed "
                    f"{i:,} / {total:,} "
                    f"assemblies "
                    f"(failed: {failed:,})"
                )


    # ========================================================
    # Final save
    # ========================================================

    save_results(
        all_rows
    )

    save_status(
        status_rows
    )


    # ========================================================
    # Summary
    # ========================================================

    result = pd.read_csv(
        OUTPUT_FILE,
        sep="\t",
    )


    status = pd.read_csv(
        STATUS_FILE,
        sep="\t",
    )


    status = status.drop_duplicates(
        subset=["genome"],
        keep="last",
    )


    n_done = (
        status["status"]
        == "done"
    ).sum()


    n_failed = (
        status["status"]
        == "failed"
    ).sum()


    print()

    print(
        "=" * 60
    )

    print(
        "DOWNLOAD COMPLETE"
    )

    print(
        "=" * 60
    )


    print(
        "Successfully queried genomes:",
        f"{n_done:,}",
    )


    print(
        "Failed genomes:",
        f"{n_failed:,}",
    )


    print(
        "Genomes with usable environmental metadata:",
        f"{result['genome'].nunique():,}",
    )


    print(
        "Environmental metadata rows:",
        f"{len(result):,}",
    )


    print()

    print(
        f"Saved:\n{OUTPUT_FILE}"
    )


    print()

    print(
        "Coverage by canonical attribute:"
    )


    coverage = (
        result
        .groupby(
            "canonical_attribute"
        )["genome"]
        .nunique()
        .sort_values(
            ascending=False
        )
    )


    print(
        coverage.to_string()
    )


    if n_failed > 0:

        print()

        print(
            f"{n_failed:,} assemblies failed. "
            "Run the script again to retry them."
        )

    else:

        print()

        print(
            "All assemblies were successfully queried."
        )


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()