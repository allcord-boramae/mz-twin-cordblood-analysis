import pandas as pd

from config import (
    MRNA_DATA2,
    MIRNA_DATA2,
    SAMPLE_METADATA,
    MRNA_INPUT_CSV,
    MIRNA_INPUT_CSV,
)

def load_sample_metadata():
    meta = pd.read_csv(SAMPLE_METADATA)
    required = {"pair_id", "T1", "T2"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"sample_metadata.csv missing columns: {sorted(missing)}")

    if len(meta) != 3:
        raise ValueError(
            f"Expected exactly 3 twin pairs, found {len(meta)}. "
            "Check data/metadata/sample_metadata.csv."
        )

    return meta

def build_sample_map(meta):
    mapping = {}
    canonical = []
    for _, row in meta.iterrows():
        pair = str(row["pair_id"]).strip()
        t1 = str(row["T1"]).strip()
        t2 = str(row["T2"]).strip()

        mapping[f"N_{t1}"] = f"{pair}-T1"
        mapping[f"N_{t2}"] = f"{pair}-T2"

        canonical.extend([f"{pair}-T1", f"{pair}-T2"])

    return mapping, canonical

def prepare_mrna(meta):
    print("=" * 70)
    print("Preparing mRNA input from Macrogen data2.xlsx")
    print("=" * 70)

    df = pd.read_excel(MRNA_DATA2)
    sample_map, canonical = build_sample_map(meta)

    required = {"ProbeID", "Gene_Symbol", *sample_map.keys()}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"mRNA data2.xlsx missing columns: {sorted(missing)}")

    cols = ["ProbeID", "Gene_Symbol", *sample_map.keys()]
    optional = [
        "Gene_ID",
        "Gene Accession",
        "Gene Description",
        "transcript_cluster_id",
        "GO Biological Process ID",
        "GO Biological Process Term",
        "GO Cellular Component ID",
        "GO Cellular Component Term",
        "GO Molecular Function ID",
        "GO Molecular Function Term",
    ]
    for col in optional:
        if col in df.columns:
            cols.append(col)

    out = df[cols].copy()
    out = out.rename(columns=sample_map)

    out["Gene_Symbol"] = (
        out["Gene_Symbol"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    out = out[out["Gene_Symbol"] != ""].copy()

    for col in canonical:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out = out.dropna(subset=canonical, how="all")

    # Preserve original row-level Macrogen feature identity.
    # For limma gene-level analysis, resolve duplicate gene symbols by
    # retaining the highest mean-expression row.
    out["mean_expr"] = out[canonical].mean(axis=1)
    out = (
        out.sort_values("mean_expr", ascending=False)
        .drop_duplicates(subset="Gene_Symbol")
        .reset_index(drop=True)
    )

    out.to_csv(MRNA_INPUT_CSV, index=False)

    print(f"Unique mRNA genes retained: {len(out):,}")
    print("Saved:", MRNA_INPUT_CSV)

def prepare_mirna(meta):
    print()
    print("=" * 70)
    print("Preparing miRNA input from Macrogen data2.xlsx")
    print("=" * 70)

    df = pd.read_excel(MIRNA_DATA2)
    sample_map, canonical = build_sample_map(meta)

    required = {"ProbeID", "mirBase", *sample_map.keys()}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"miRNA data2.xlsx missing columns: {sorted(missing)}")

    cols = ["ProbeID", "mirBase", *sample_map.keys()]
    optional = [
        "Probe Set Name",
        "Accession",
        "Transcript ID(Array Design)",
        "GeneChip Array",
        "Species Scientific Name",
        "Annotation Date",
    ]
    for col in optional:
        if col in df.columns:
            cols.append(col)

    out = df[cols].copy()
    out = out.rename(columns=sample_map)

    out["mirBase"] = (
        out["mirBase"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    out = out[out["mirBase"].str.startswith("hsa-")].copy()

    for col in canonical:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out = out.dropna(subset=canonical, how="all")
    out["mean_expr"] = out[canonical].mean(axis=1)

    out = (
        out.sort_values("mean_expr", ascending=False)
        .drop_duplicates(subset="mirBase")
        .reset_index(drop=True)
    )

    out.to_csv(MIRNA_INPUT_CSV, index=False)

    print(f"Unique human miRNAs retained: {len(out):,}")
    print("Saved:", MIRNA_INPUT_CSV)

def main():
    for p in [MRNA_DATA2, MIRNA_DATA2, SAMPLE_METADATA]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input: {p}")

    meta = load_sample_metadata()

    print("Twin pairs used for analysis:")
    print(meta.to_string(index=False))
    print()

    prepare_mrna(meta)
    prepare_mirna(meta)

    print()
    print("Input preparation completed successfully.")

if __name__ == "__main__":
    main()
