#!/bin/bash
set -euo pipefail

echo "[1/4] Prepare analysis inputs from Macrogen data2.xlsx"
python3 src/01_prepare_inputs.py

echo "[2/4] Paired limma + ordinary paired t-test + BH-FDR"
Rscript src/02_paired_limma_DE.R   data/processed/mRNA_analysis_input.csv   data/processed/miRNA_analysis_input.csv   results

echo "[3/4] Original-vs-revised audit"
python3 src/R01_original_vs_revised_DE.py

echo "[4/4] Verify manuscript dual-test candidate sets"
python3 src/verify_dual_test_outputs.py results

echo
echo "DONE"
echo "Main results: results/"
echo "Revision audit: results/revision/"
