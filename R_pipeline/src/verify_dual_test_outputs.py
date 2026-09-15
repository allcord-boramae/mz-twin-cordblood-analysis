"""Verify that the canonical R outputs implement the manuscript dual-test rule."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


FC_THRESHOLD = 1.5
P_THRESHOLD = 0.05

EXPECTED_DEMS = {
    "hsa-miR-6779-5p",
    "hsa-miR-1292-5p",
    "hsa-miR-6891-5p",
    "hsa-miR-4253",
    "hsa-miR-128-3p",
    "hsa-miR-148b-3p",
    "hsa-miR-3679-5p",
    "hsa-miR-328-3p",
}

EXPECTED_DEGS = {
    "IGLV1-40",
    "LOC100128751",
    "PLEKHM1P",
    "FAM91A1",
    "LOC100131174",
    "UBE2Q2L",
    "SLED1",
    "TOMM40L",
    "PHC1",
    "IGKV2-29",
    "SH2D1A",
    "MIR339",
    "RABGGTB",
    "DNTTIP2",
    "CISD1",
    "LOC440602",
    "SLC52A3",
    "R3HDM1",
    "SGOL1",
    "LOC100128653",
    "PRKX",
    "OR2L3",
    "LOC100996579",
}


def as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().eq("true")


def verify_dataset(
    results_dir: Path,
    *,
    full_name: str,
    candidate_name: str,
    moderated_only_name: str,
    feature_column: str,
    expected_candidates: set[str],
    expected_moderated_only_count: int,
) -> None:
    full = pd.read_csv(results_dir / full_name)
    candidates = pd.read_csv(results_dir / candidate_name)
    moderated_only = pd.read_csv(results_dir / moderated_only_name)

    required = {
        feature_column,
        "logFC",
        "P.Value",
        "P.ordinary_t",
        "consistent_direction",
        "moderated_only_candidate",
        "exploratory_candidate",
    }
    missing = required - set(full.columns)
    if missing:
        raise AssertionError(f"{full_name} missing columns: {sorted(missing)}")

    consistent = as_bool(full["consistent_direction"])
    dual_mask = (
        (2 ** full["logFC"].abs() >= FC_THRESHOLD)
        & (full["P.Value"] < P_THRESHOLD)
        & (full["P.ordinary_t"] < P_THRESHOLD)
        & consistent
    )
    moderated_only_mask = (
        (2 ** full["logFC"].abs() >= FC_THRESHOLD)
        & (full["P.Value"] < P_THRESHOLD)
        & consistent
    )

    if not as_bool(full["exploratory_candidate"]).equals(dual_mask):
        raise AssertionError(f"{full_name}: exploratory_candidate does not match the dual-test rule")
    if not as_bool(full["moderated_only_candidate"]).equals(moderated_only_mask):
        raise AssertionError(
            f"{full_name}: moderated_only_candidate does not match the moderated-only rule"
        )

    calculated_candidates = set(full.loc[dual_mask, feature_column])
    written_candidates = set(candidates[feature_column])
    if calculated_candidates != written_candidates:
        raise AssertionError(f"{candidate_name} does not match the full-results flags")
    if written_candidates != expected_candidates:
        missing_expected = sorted(expected_candidates - written_candidates)
        unexpected = sorted(written_candidates - expected_candidates)
        raise AssertionError(
            f"{candidate_name} manuscript set mismatch; "
            f"missing={missing_expected}, unexpected={unexpected}"
        )

    calculated_moderated_only = set(full.loc[moderated_only_mask, feature_column])
    written_moderated_only = set(moderated_only[feature_column])
    if calculated_moderated_only != written_moderated_only:
        raise AssertionError(f"{moderated_only_name} does not match the full-results flags")
    if len(written_moderated_only) != expected_moderated_only_count:
        raise AssertionError(
            f"{moderated_only_name}: expected {expected_moderated_only_count}, "
            f"found {len(written_moderated_only)}"
        )

    print(
        f"{feature_column}: {len(written_candidates)} dual-test candidates; "
        f"{len(written_moderated_only)} moderated-test-only candidates"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    args = parser.parse_args()

    verify_dataset(
        args.results_dir,
        full_name="miRNA_limma_full_results.csv",
        candidate_name="DEM_exploratory_candidates.csv",
        moderated_only_name="DEM_moderated_only_candidates.csv",
        feature_column="miRNA",
        expected_candidates=EXPECTED_DEMS,
        expected_moderated_only_count=30,
    )
    verify_dataset(
        args.results_dir,
        full_name="mRNA_limma_full_results.csv",
        candidate_name="DEG_exploratory_candidates.csv",
        moderated_only_name="DEG_moderated_only_candidates.csv",
        feature_column="Gene",
        expected_candidates=EXPECTED_DEGS,
        expected_moderated_only_count=114,
    )
    print("Dual-test output verification passed.")


if __name__ == "__main__":
    main()
