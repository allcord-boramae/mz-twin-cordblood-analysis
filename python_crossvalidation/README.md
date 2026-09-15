# Integrated miRNA–mRNA expression profiling in cord blood from monozygotic twins — analysis code

> **Note (2026-09-13).** The canonical statistics reported in the manuscript
> are produced by the R/Bioconductor limma pipeline (limma v3.68.5, R 4.6.1;
> see the accompanying `R_pipeline` directory and its `R_sessionInfo_limma.txt`).
> The Python pipeline in this repository is an independent re-implementation
> used for cross-validation; it reproduces the candidate lists and all
> conclusions, with small numerical differences in the moderated statistics
> attributable to the empirical-Bayes prior estimation.

Analysis code accompanying:

> Roh EY, Oh I, Shin S, Yoon JH, Kim BJ. *Integrated miRNA–mRNA expression
> profiling in cord blood from three pairs of monozygotic twins.* PLOS ONE
> (PONE-D-26-14953).

All analyses run directly on the RMA-normalized expression matrices deposited
in NCBI GEO, so every result in the manuscript can be reproduced from the
public data alone.

## Data

Download the processed data files from GEO and place them in a local
directory (referred to as `<GEO_DIR>` below):

| GEO accession | Platform | File |
|---|---|---|
| [GSE326366](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326366) | Affymetrix Human Gene 2.0 ST (GPL16686) | `processed_mRNA_expression_matrix.txt` |
| [GSE326361](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326361) | Affymetrix miRNA 4.0 (GPL19117) | `processed_miRNA_expression_matrix.txt` |

Samples: three monozygotic twin pairs (MZ1–MZ3); T1 = first-born,
T2 = second-born. Both matrices contain the same six samples
(`MZ1-T1` … `MZ3-T2`).

The GEO matrices alone do not carry probe/GO annotation, so the repository
ships the two array annotation files the pipeline needs (both derived from
the public Affymetrix annotation releases):

- `annotation/miRNA40_probe_annotation.tsv` — Affymetrix miRNA 4.0 probe ID
  → miRBase v20 miRNA name mapping (used to restrict the miRNA analysis to
  human `hsa-` miRNAs; the mRNA matrix already contains gene symbols).
- `annotation/HuGene20ST_GO_annotation.tsv.gz` — Human Gene 2.0 ST gene
  symbol → GO term (BP/CC/MF) mapping, used by the GO enrichment script.

Provenance: both files were extracted verbatim from the official Affymetrix
(Thermo Fisher) NetAffx array annotation as distributed with the original
array analysis deliverable, and are frozen here so that the pipeline uses
exactly the annotation version underlying the published results. The same
annotation is independently available from the Thermo Fisher annotation
downloads and from Bioconductor (`hugene20sttranscriptcluster.db`,
`pd.mirna.4.0`), so the files can be verified against those public sources.
Raw .CEL files are also deposited in GEO, allowing full re-analysis from
scratch (RMA normalization onward) if desired.

## Requirements

Python ≥ 3.10 with `numpy`, `pandas`, `matplotlib`, `Pillow`, and `scipy` (≥ 1.11 for
`scipy.stats.false_discovery_control`).

```
pip install numpy pandas scipy matplotlib Pillow
```

## Scripts

### 1. Differential expression (`scripts/01_differential_expression.py`)

Paired T2-vs-T1 analysis of within-pair log2 differences using a
limma-equivalent empirical Bayes moderated t-test (Smyth 2004; the
`fit_f_dist`, `trigamma_inverse` and `squeeze_var` functions are direct
ports of `limma::fitFDist`, `limma::trigammaInverse` and
`limma::squeezeVar`), with Benjamini–Hochberg FDR correction over all tested
features. Candidate DEMs/DEGs are selected by a high-stringency dual-test
criterion: |FC| ≥ 1.5, the same direction of change in all three pairs, and
P < 0.05 under both the moderated t-test and the ordinary paired t-test
(8 DEMs, 23 DEGs). The wider moderated-test-only lists and full per-feature
results with FDR-adjusted P values are written alongside. Also runs the
female-pair (MZ1, MZ2) sensitivity analysis.

```
python scripts/01_differential_expression.py \
    --geo-dir <GEO_DIR> \
    --annotation annotation/miRNA40_probe_annotation.tsv \
    --out results
```

Outputs: `miRNA_moderated_all_results.csv`, `Gene_moderated_all_results.csv`,
`DEM_candidates.csv`, `DEG_candidates.csv`.

### 2. Cell-type deconvolution (`scripts/02_deconvolution.py`)

NNLS (CIBERSORT-like) estimation of 11 cord blood cell-type proportions from
the mRNA matrix, followed by within-pair comparisons using the exact paired
Wilcoxon signed-rank test. With n = 3 pairs the minimum attainable exact
two-sided P value is 0.25, so these comparisons are descriptive.

```
python scripts/02_deconvolution.py --geo-dir <GEO_DIR> --out results
```

Outputs: `deconvolution_NNLS_proportions.csv`,
`deconvolution_within_pair_stats.csv`.

### 3. Figures (`scripts/03_figures.py`, `scripts/04_supplementary_figures.py`)

Regenerate the main figures (Fig 1–4) and the supplementary figures (S2–S3
Fig) as PLOS-compliant LZW TIFFs.

### 4. GO enrichment (`scripts/05_go_enrichment.py`)

One-sided Fisher's exact test for the 23 candidate DEGs against the
29,383-gene array background. As in standard enrichment tools, only GO terms
annotated to at least one candidate DEG are testable (m = 126); BH FDR
correction is applied across these tests, and full per-term results are
written to `results/GO_enrichment_full_results.csv`.

```
python scripts/05_go_enrichment.py --geo-dir <GEO_DIR> \
    --annotation annotation/HuGene20ST_GO_annotation.tsv.gz \
    --deg results/DEG_candidates.csv --out results
```

## License

MIT (code). The expression data are available from NCBI GEO under the
accessions above.
