#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript src/02_paired_limma_DE.R <mRNA_input.csv> <miRNA_input.csv> <results_dir>")
}

mrna_file <- args[1]
mirna_file <- args[2]
results_dir <- args[3]

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(
  "MZ1-T1", "MZ1-T2",
  "MZ2-T1", "MZ2-T2",
  "MZ3-T1", "MZ3-T2"
)

fc_threshold <- log2(1.5)

read_expression <- function(path, feature_col) {
  df <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required <- c(feature_col, samples)
  missing <- setdiff(required, colnames(df))

  if (length(missing) > 0) {
    stop(
      paste(
        basename(path),
        "missing columns:",
        paste(missing, collapse = ", ")
      )
    )
  }

  df <- df[
    !is.na(df[[feature_col]]) &
    df[[feature_col]] != "",
  ]

  df <- df[
    !duplicated(df[[feature_col]]),
  ]

  rownames(df) <- make.unique(
    as.character(df[[feature_col]])
  )

  as.matrix(
    df[, samples, drop = FALSE]
  )
}

run_paired_limma <- function(expr, feature_name) {

  # One independent observation per twin pair:
  # within-pair T2 - T1 expression difference.
  diffs <- cbind(
    MZ1 = expr[, "MZ1-T2"] - expr[, "MZ1-T1"],
    MZ2 = expr[, "MZ2-T2"] - expr[, "MZ2-T1"],
    MZ3 = expr[, "MZ3-T2"] - expr[, "MZ3-T1"]
  )

  design <- matrix(
    1,
    nrow = 3,
    ncol = 1
  )

  colnames(design) <- "T2_minus_T1"

  fit <- lmFit(
    diffs,
    design
  )

  fit <- eBayes(
    fit,
    robust = TRUE
  )

  tt <- topTable(
    fit,
    coef = "T2_minus_T1",
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )

  # Ordinary paired t-test, expressed as a one-sample t-test on the three
  # within-pair differences. Zero-variance rows are assigned P = 1, matching
  # the conservative rule used by the independent Python implementation.
  ordinary_p <- apply(
    diffs,
    1,
    function(x) {
      s <- sd(x)

      if (!is.finite(s) || s <= 0) {
        return(1)
      }

      t_stat <- mean(x) / (s / sqrt(length(x)))
      2 * pt(-abs(t_stat), df = length(x) - 1)
    }
  )

  result <- data.frame(
    feature = rownames(diffs),
    logFC = rowMeans(diffs),
    FC_abs = 2 ^ abs(rowMeans(diffs)),
    diff_MZ1 = diffs[, "MZ1"],
    diff_MZ2 = diffs[, "MZ2"],
    diff_MZ3 = diffs[, "MZ3"],
    t_moderated = tt$t,
    P.Value = tt$P.Value,
    adj.P.Val = tt$adj.P.Val,
    P.ordinary_t = ordinary_p,
    B = tt$B,
    stringsAsFactors = FALSE
  )

  result$consistent_direction <- apply(
    diffs,
    1,
    function(x) all(x > 0) || all(x < 0)
  )

  result$FDR_significant <- (
    result$adj.P.Val < 0.05
  )

  # Wider moderated-test-only list retained for revision audits.
  result$moderated_only_candidate <- (
    abs(result$logFC) >= fc_threshold &
    result$P.Value < 0.05 &
    result$consistent_direction
  )

  # Manuscript exploratory dual-test rule.
  result$exploratory_candidate <- (
    abs(result$logFC) >= fc_threshold &
    result$P.Value < 0.05 &
    result$P.ordinary_t < 0.05 &
    result$consistent_direction
  )

  colnames(result)[1] <- feature_name

  result[
    order(result$P.Value),
  ]
}

mrna <- read_expression(
  mrna_file,
  "Gene_Symbol"
)

mirna <- read_expression(
  mirna_file,
  "mirBase"
)

mrna_res <- run_paired_limma(
  mrna,
  "Gene"
)

mirna_res <- run_paired_limma(
  mirna,
  "miRNA"
)

write.csv(
  mrna_res,
  file.path(results_dir, "mRNA_limma_full_results.csv"),
  row.names = FALSE
)

write.csv(
  mirna_res,
  file.path(results_dir, "miRNA_limma_full_results.csv"),
  row.names = FALSE
)

write.csv(
  mrna_res[mrna_res$FDR_significant, ],
  file.path(results_dir, "mRNA_FDR_significant.csv"),
  row.names = FALSE
)

write.csv(
  mirna_res[mirna_res$FDR_significant, ],
  file.path(results_dir, "miRNA_FDR_significant.csv"),
  row.names = FALSE
)

write.csv(
  mrna_res[mrna_res$exploratory_candidate, ],
  file.path(results_dir, "DEG_exploratory_candidates.csv"),
  row.names = FALSE
)

write.csv(
  mirna_res[mirna_res$exploratory_candidate, ],
  file.path(results_dir, "DEM_exploratory_candidates.csv"),
  row.names = FALSE
)

write.csv(
  mrna_res[mrna_res$moderated_only_candidate, ],
  file.path(results_dir, "DEG_moderated_only_candidates.csv"),
  row.names = FALSE
)

write.csv(
  mirna_res[mirna_res$moderated_only_candidate, ],
  file.path(results_dir, "DEM_moderated_only_candidates.csv"),
  row.names = FALSE
)

summary <- data.frame(
  dataset = c("mRNA", "miRNA"),
  n_tested = c(
    nrow(mrna_res),
    nrow(mirna_res)
  ),
  n_FDR_lt_0.05 = c(
    sum(mrna_res$FDR_significant),
    sum(mirna_res$FDR_significant)
  ),
  n_exploratory_candidates = c(
    sum(mrna_res$exploratory_candidate),
    sum(mirna_res$exploratory_candidate)
  ),
  n_moderated_only_candidates = c(
    sum(mrna_res$moderated_only_candidate),
    sum(mirna_res$moderated_only_candidate)
  ),
  min_FDR = c(
    min(mrna_res$adj.P.Val, na.rm = TRUE),
    min(mirna_res$adj.P.Val, na.rm = TRUE)
  )
)

write.csv(
  summary,
  file.path(results_dir, "DE_summary.csv"),
  row.names = FALSE
)

print(summary)

sink(
  file.path(results_dir, "R_sessionInfo_limma.txt")
)

print(sessionInfo())

sink()
