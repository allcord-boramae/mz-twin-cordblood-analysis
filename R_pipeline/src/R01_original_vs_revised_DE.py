import pandas as pd

from config import (
    MRNA_LIMMA_FULL,
    MIRNA_LIMMA_FULL,
    REVISION_RESULTS_DIR,
)

ORIGINAL_DEMS = [
    "hsa-miR-6779-5p",
    "hsa-miR-1292-5p",
    "hsa-miR-6891-5p",
    "hsa-miR-4253",
    "hsa-miR-128-3p",
    "hsa-miR-148b-3p",
    "hsa-miR-3679-5p",
    "hsa-miR-328-3p",
]

ORIGINAL_DEGS = [
    "CISD1","TOMM40L","PRKX","PHC1","SGOL1","TNKS2","ARMC8","MGAT5B",
    "ZNF33A","ZNF585B","NCOA1","DLC1","ZNF212","CPTP","TIMM17B","PPM1G",
    "USP9Y","NUPL2","MED13L","ZNF544","ERCC8","ANKRD36","PDZD8",
]

def status_label(row):
    if bool(row["FDR_significant"]):
        return "FDR-significant"
    if bool(row["exploratory_candidate"]):
        return "Exploratory only"
    return "Not retained"

def main():
    for p in [MRNA_LIMMA_FULL, MIRNA_LIMMA_FULL]:
        if not p.exists():
            raise FileNotFoundError(p)

    mi = pd.read_csv(MIRNA_LIMMA_FULL)
    mr = pd.read_csv(MRNA_LIMMA_FULL)

    dem = mi[mi["miRNA"].isin(ORIGINAL_DEMS)].copy()
    deg = mr[mr["Gene"].isin(ORIGINAL_DEGS)].copy()

    dem["Revised status"] = dem.apply(status_label, axis=1)
    deg["Revised status"] = deg.apply(status_label, axis=1)

    summary = pd.DataFrame([
        {
            "Original set": "8 DEMs",
            "Found": len(dem),
            "Still exploratory": int(dem["exploratory_candidate"].sum()),
            "BH-FDR significant": int(dem["FDR_significant"].sum()),
        },
        {
            "Original set": "23 DEGs",
            "Found": len(deg),
            "Still exploratory": int(deg["exploratory_candidate"].sum()),
            "BH-FDR significant": int(deg["FDR_significant"].sum()),
        },
    ])

    REVISION_RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    dem.to_csv(
        REVISION_RESULTS_DIR / "R01_original_8_DEM_vs_revised_limma.csv",
        index=False
    )

    deg.to_csv(
        REVISION_RESULTS_DIR / "R01_original_23_DEG_vs_revised_limma.csv",
        index=False
    )

    summary.to_csv(
        REVISION_RESULTS_DIR / "R01_revision_summary.csv",
        index=False
    )

    with pd.ExcelWriter(
        REVISION_RESULTS_DIR / "R01_original_vs_revised_DE.xlsx"
    ) as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        dem.to_excel(writer, sheet_name="Original 8 DEMs", index=False)
        deg.to_excel(writer, sheet_name="Original 23 DEGs", index=False)

    print(summary.to_string(index=False))
    print()
    print("Saved revision audit to:", REVISION_RESULTS_DIR)

if __name__ == "__main__":
    main()
