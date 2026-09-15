# mz-twin-cordblood-analysis

Analysis code accompanying:

> Roh EY*, Oh I*, Park J*, Shin S, Yoon JH, Kim BJ. *Integrated miRNA-mRNA expression
> profiling in cord blood from three pairs of monozygotic twins.* PLOS ONE
> (PONE-D-26-14953). (*co-first authors)

All analyses run on the RMA-normalized expression matrices deposited in NCBI GEO
(GSE326361 miRNA, GSE326366 mRNA; released upon publication, reviewer tokens during
review), so every result in the manuscript can be reproduced from the deposited data
plus the annotation files in this repository.

## Layout

| Folder | Contents |
|---|---|
| `R_pipeline/` | **Canonical statistics** reported in the manuscript: paired differential expression with Bioconductor limma (v3.68.5, R 4.6.1; `results/R_sessionInfo_limma.txt`), ordinary paired t-tests, BH FDR correction, and dual-test candidate selection (8 DEMs, 23 DEGs). The wider moderated-test-only sets (30 DEMs, 114 DEGs) are written separately for revision audits. Run via `run_all.sh` or the numbered notebooks. |
| `python_crossvalidation/` | Independent Python re-implementation used for cross-validation (reproduces the candidate lists and all conclusions), plus the cell-type deconvolution (NNLS, marker z-score, ssGSEA; exact paired Wilcoxon), GO enrichment, and figure-generation scripts. See its own README for usage. |
| `annotation/` | Array annotation required by the pipelines: Affymetrix miRNA 4.0 probe → miRBase v20 mapping and Human Gene 2.0 ST gene → GO mapping (both extracted verbatim from the official Affymetrix NetAffx annotation; independently verifiable against Thermo Fisher / Bioconductor releases). |
| `data_reference/TargetScan8.0_exports/` | TargetScan 8.0 per-miRNA predicted-target exports (2026-03-22) used to build and verify the miRNA-mRNA cross-validation table (Table 4), including non-conserved miRNA families absent from the TargetScan summary files. |

**Note on the `R_pipeline` revision notebooks (R01–R03) and tracked integration outputs.**
The reviewer-response R02–R03 notebooks explicitly use the *wider audit lists*
(30 DEMs × 114 DEGs, moderated-test-only). The already tracked integration and GO
outputs inside `R_pipeline` (`results/miRNA_mRNA_*`, `results/revision/`, and
`results/GO_enrichment_full_results.csv`) were also generated from those wider
lists and, for TargetScan, from the conserved-family summary file only; they
therefore differ by design from the manuscript's Table 4 and S3 Table. The standard
R notebooks now consume the dual-test candidate files. The numbers reported in the manuscript use
the dual-test candidate lists (8 DEMs × 23 DEGs) and were verified directly
against the database source files — the per-miRNA TargetScan 8.0 exports in
`data_reference/` (which include non-conserved families), the miRDB v6.0
prediction set (matched via RefSeq accessions), and miRTarBase v9.0 MTI
records; the GO enrichment in the manuscript tests the 23 candidate DEGs
against the full 29,383-gene array background (`python_crossvalidation/scripts/05_go_enrichment.py`).
The R/limma statistics (`R_pipeline/results/*_limma_full_results.csv`) are the
canonical source of the moderated P values reported in the manuscript. The
Python implementation independently reproduces the candidate sets and
conclusions, but its empirical-Bayes P values can differ slightly.

Expression data are **not** stored here — download the processed matrices from GEO
(see `python_crossvalidation/README.md` for step-by-step instructions).
The full miRDB v6.0 and miRTarBase v9.0 files are not redistributed here
(size/licensing); they are freely available from http://mirdb.org and
https://mirtarbase.cuhk.edu.cn.

## License

MIT (code). Array annotation: Affymetrix/Thermo Fisher public annotation releases.
TargetScan exports: TargetScan v8.0 (CC BY-NC 4.0, Whitehead Institute).
