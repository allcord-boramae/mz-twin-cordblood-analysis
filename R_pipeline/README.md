# MCDA Revision Final

This project follows the original study workflow confirmed by the first author:

1. Raw Affymetrix CEL files were generated.
2. Macrogen performed RMA normalization and platform annotation.
3. Macrogen supplied `data2.xlsx`.
4. The study team performed downstream statistical analysis from those processed files.

Therefore, the **revision statistical analysis starts from the original Macrogen `data2.xlsx` files**.
Raw CEL files are retained separately for raw-data submission / GEO availability.

## Required files

### Statistical-analysis inputs
Copy:

- mRNA `data2.xlsx` → `data/processed/mRNA/data2.xlsx`
- miRNA `data2.xlsx` → `data/processed/miRNA/data2.xlsx`

### Raw-data submission files
Put raw CEL files under:

- `data/raw/mRNA_CEL/`
- `data/raw/miRNA_CEL/`

The CEL files are not re-normalized in this revision pipeline.

## Sample pairs used

The current project is configured according to the manuscript analysis:

- MZ1: 2041C / 2042C
- MZ2: 3640C / 3641C
- MZ3: 12156C / 12157C

These are defined in:

`data/metadata/sample_metadata.csv`

Do not hard-code sample IDs elsewhere.

## Statistical revision

Original approach:
- paired t-test
- FC >= 1.5
- nominal P < 0.05
- same direction in all 3 pairs
- no DE-level multiple-testing correction

Revised approach:
- Bioconductor `limma`
- within-pair T2-T1 differences for the 3 independent twin pairs
- empirical-Bayes moderated statistics
- ordinary paired t-test on the same three within-pair differences
- Benjamini-Hochberg FDR across all tested features

Confirmatory:
- BH-FDR < 0.05

Exploratory:
- |log2FC| >= log2(1.5)
- moderated P < 0.05
- ordinary paired t-test P < 0.05
- same direction in all 3 pairs

This dual-test rule yields the manuscript candidate sets: 8 DEMs and 23 DEGs.
The wider moderated-test-only sets (30 DEMs and 114 DEGs) are retained as
separate `*_moderated_only_candidates.csv` audit outputs.

The exploratory rule does not replace multiple-testing correction.

## Results organization

Main analysis:
`results/`

- `DEM_exploratory_candidates.csv`: 8 dual-test DEMs
- `DEG_exploratory_candidates.csv`: 23 dual-test DEGs
- `DEM_moderated_only_candidates.csv`: 30 moderated-test-only DEMs
- `DEG_moderated_only_candidates.csv`: 114 moderated-test-only DEGs

Reviewer-response / audit outputs:
`results/revision/`

## Run

From project root:

```bash
python3 -m pip install -r requirements.txt
bash run_all.sh
```

Or run the notebooks in order:

- `00_setup.ipynb`
- `01_input_qc.ipynb`
- `02_differential_expression.ipynb`
- `revision/R01_original_vs_revised_DE.ipynb`

Downstream notebooks are intentionally left controlled until the revised DE result is reconfirmed.
