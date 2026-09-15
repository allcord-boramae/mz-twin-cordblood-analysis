from pathlib import Path

GLOBAL_SEED = 20260829

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
METADATA_DIR = DATA_DIR / "metadata"

MRNA_CEL_DIR = RAW_DIR / "mRNA_CEL"
MIRNA_CEL_DIR = RAW_DIR / "miRNA_CEL"

MRNA_DATA2 = PROCESSED_DIR / "mRNA" / "data2.xlsx"
MIRNA_DATA2 = PROCESSED_DIR / "miRNA" / "data2.xlsx"

SAMPLE_METADATA = METADATA_DIR / "sample_metadata.csv"

RESULTS_DIR = PROJECT_ROOT / "results"
REVISION_RESULTS_DIR = RESULTS_DIR / "revision"
FIGURES_DIR = PROJECT_ROOT / "figures"
REVISION_FIGURES_DIR = FIGURES_DIR / "revision"
LOG_DIR = PROJECT_ROOT / "logs"

for d in [
    MRNA_CEL_DIR, MIRNA_CEL_DIR,
    PROCESSED_DIR, METADATA_DIR,
    RESULTS_DIR, REVISION_RESULTS_DIR,
    FIGURES_DIR, REVISION_FIGURES_DIR,
    LOG_DIR
]:
    d.mkdir(parents=True, exist_ok=True)

# Standardized analysis outputs
MRNA_INPUT_CSV = PROCESSED_DIR / "mRNA_analysis_input.csv"
MIRNA_INPUT_CSV = PROCESSED_DIR / "miRNA_analysis_input.csv"

# Full DE results
MRNA_LIMMA_FULL = RESULTS_DIR / "mRNA_limma_full_results.csv"
MIRNA_LIMMA_FULL = RESULTS_DIR / "miRNA_limma_full_results.csv"

MRNA_FDR = RESULTS_DIR / "mRNA_FDR_significant.csv"
MIRNA_FDR = RESULTS_DIR / "miRNA_FDR_significant.csv"

DEG_EXPLORATORY = RESULTS_DIR / "DEG_exploratory_candidates.csv"
DEM_EXPLORATORY = RESULTS_DIR / "DEM_exploratory_candidates.csv"
DEG_MODERATED_ONLY = RESULTS_DIR / "DEG_moderated_only_candidates.csv"
DEM_MODERATED_ONLY = RESULTS_DIR / "DEM_moderated_only_candidates.csv"

DE_SUMMARY = RESULTS_DIR / "DE_summary.csv"
