"""
Differential expression analysis of MZ twin cord blood (T2 vs T1, paired).

Input:  RMA-normalized expression matrices exactly as deposited in NCBI GEO
        (GSE326361 miRNA, GSE326366 mRNA):
          processed_miRNA_expression_matrix.txt  (ID_REF + MZ1-T1 ... MZ3-T2)
          processed_mRNA_expression_matrix.txt   (ID_REF, Gene_Symbol + samples)
        plus annotation/miRNA40_probe_annotation.tsv (Affymetrix GeneChip
        miRNA 4.0 probe -> miRBase v20 name mapping).

Statistics: within-pair log2 differences (T2 - T1) are fitted with a
limma-equivalent empirical Bayes moderated paired t-test (Smyth 2004;
the fit_f_dist / trigamma_inverse / squeeze_var functions are direct ports
of limma::fitFDist, limma::trigammaInverse and limma::squeezeVar), followed
by Benjamini-Hochberg FDR correction across all tested features.

Selection of candidate DEMs/DEGs (exploratory, as stated in the manuscript):
a high-stringency dual-test criterion requiring |logFC| >= log2(1.5),
identical direction of change in all three twin pairs, and P < 0.05 under
BOTH the empirical Bayes moderated t-test and the ordinary paired t-test.
The full moderated-test results (all features, with BH FDR-adjusted P
values) and the wider single-test candidate lists are also written out as
supplementary outputs.

Also runs the sensitivity analysis restricted to the two female pairs
(MZ1, MZ2) and reports overlap with the three-pair candidate lists.

Usage:
  python 01_differential_expression.py --geo-dir <dir with the two matrices>
      --annotation ../annotation/miRNA40_probe_annotation.tsv --out ../results

Requires: python >= 3.10, numpy, pandas, scipy.
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.special import digamma, polygamma
from scipy import stats

T1_COLS = ['MZ1-T1', 'MZ2-T1', 'MZ3-T1']
T2_COLS = ['MZ1-T2', 'MZ2-T2', 'MZ3-T2']
FEMALE_IDX = [0, 1]          # MZ1, MZ2 are the female pairs
FC_THR = np.log2(1.5)
P_THR = 0.05


# ---------------------------------------------------------------------------
# limma ports (Smyth 2004, Stat Appl Genet Mol Biol 3:Article3)
# ---------------------------------------------------------------------------
def trigamma(x):
    return polygamma(1, x)


def trigamma_inverse(x):
    """Port of limma::trigammaInverse (Newton iteration)."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = x > 1e7
    out[small] = 1.0 / np.sqrt(x[small])
    big = x < 1e-6
    out[big] = 1.0 / x[big]
    mid = ~(small | big)
    y = 0.5 + 1.0 / x[mid]
    for _ in range(50):
        tri = trigamma(y)
        dif = tri * (1 - tri / x[mid]) / polygamma(2, y)
        y = y + dif
        if np.max(-dif / y) < 1e-8:
            break
    out[mid] = y
    return out


def fit_f_dist(s2, df1):
    """Port of limma::fitFDist. Returns (prior variance s0^2, prior df d0)."""
    s2 = np.asarray(s2, dtype=float)
    x = s2[np.isfinite(s2) & (s2 > 0)]
    z = np.log(x)
    e = z - digamma(df1 / 2.0) + np.log(df1 / 2.0)
    emean = e.mean()
    evar = np.sum((e - emean) ** 2) / (e.size - 1) - trigamma(df1 / 2.0)
    if evar > 0:
        d0 = 2 * trigamma_inverse(np.array([evar]))[0]
        s20 = np.exp(emean + digamma(d0 / 2.0) - np.log(d0 / 2.0))
    else:
        d0 = np.inf
        s20 = np.exp(emean)
    return s20, d0


def squeeze_var(s2, df):
    """Port of limma::squeezeVar."""
    s20, d0 = fit_f_dist(s2, df)
    if np.isinf(d0):
        return np.full_like(s2, s20), s20, d0
    return (d0 * s20 + df * s2) / (d0 + df), s20, d0


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def moderated_paired_analysis(expr, feature_name):
    """Moderated paired T2-vs-T1 analysis of a features x 6-sample matrix."""
    diffs = expr[T2_COLS].to_numpy() - expr[T1_COLS].to_numpy()
    n = diffs.shape[1]
    df_resid = n - 1
    logfc = diffs.mean(axis=1)
    s2 = diffs.var(axis=1, ddof=1)

    s2_post, s20, d0 = squeeze_var(s2, df_resid)
    df_total = d0 + df_resid
    t_mod = logfc / (np.sqrt(s2_post) / np.sqrt(n))
    if np.isinf(df_total):
        p_mod = 2 * stats.norm.sf(np.abs(t_mod))
    else:
        p_mod = 2 * stats.t.sf(np.abs(t_mod), df=df_total)
    p_adj = stats.false_discovery_control(p_mod, method='bh')

    # ordinary paired t-test for the dual-test criterion
    with np.errstate(divide='ignore', invalid='ignore'):
        t_ord = logfc / np.sqrt(s2 / n)
    p_ord = np.where(s2 > 0, 2 * stats.t.sf(np.abs(t_ord), df=df_resid), 1.0)

    signs = np.sign(diffs)
    consistent = (np.abs(signs.sum(axis=1)) == n)

    res = pd.DataFrame({
        feature_name: expr.index,
        'logFC': logfc,
        'diff_MZ1': diffs[:, 0],
        'diff_MZ2': diffs[:, 1],
        'diff_MZ3': diffs[:, 2],
        't_moderated': t_mod,
        'P.Value': p_mod,
        'adj.P.Val': p_adj,
        'P.ordinary_t': p_ord,
        'consistent_direction': consistent,
    }).sort_values('P.Value').reset_index(drop=True)
    return res, {'prior_df_d0': d0, 'prior_var_s02': s20}


def select_candidates(res):
    """High-stringency dual-test criterion (see module docstring)."""
    return res[(res['P.Value'] < P_THR)
               & (res['P.ordinary_t'] < P_THR)
               & (res['logFC'].abs() >= FC_THR)
               & res['consistent_direction']].copy()


def select_moderated_only(res):
    """Wider list: moderated-t criterion alone (supplementary output)."""
    return res[(res['P.Value'] < P_THR)
               & (res['logFC'].abs() >= FC_THR)
               & res['consistent_direction']].copy()


def female_pair_sensitivity(expr):
    """Two-female-pair analysis: |mean logFC| >= threshold, same direction."""
    d = expr[[T2_COLS[i] for i in FEMALE_IDX]].to_numpy() \
        - expr[[T1_COLS[i] for i in FEMALE_IDX]].to_numpy()
    logfc = d.mean(axis=1)
    same_dir = (np.sign(d[:, 0]) == np.sign(d[:, 1])) & (d != 0).all(axis=1)
    keep = (np.abs(logfc) >= FC_THR) & same_dir
    return set(expr.index[keep])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geo-dir', required=True)
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # miRNA: annotate probes, keep human miRNAs, deduplicate
    mi = pd.read_csv(os.path.join(args.geo_dir, 'processed_miRNA_expression_matrix.txt'), sep='\t')
    ann = pd.read_csv(args.annotation, sep='\t')
    mi = mi.merge(ann, left_on='ID_REF', right_on='ProbeID', how='left')
    mi['miRNA_name'] = mi['miRNA_name'].fillna(mi['ID_REF'].astype(str))
    mi = mi[mi['miRNA_name'].str.startswith('hsa-')]
    expr_mi = mi.set_index('miRNA_name')[T1_COLS + T2_COLS]
    expr_mi = expr_mi[~expr_mi.index.duplicated(keep='first')]
    print(f"miRNA: {len(expr_mi)} unique human (hsa-) features")

    # mRNA: keep annotated genes, collapse duplicates to highest-mean probe
    mr = pd.read_csv(os.path.join(args.geo_dir, 'processed_mRNA_expression_matrix.txt'), sep='\t')
    mr = mr[mr['Gene_Symbol'].notna()].copy()
    mr['mean_expr'] = mr[T1_COLS + T2_COLS].mean(axis=1)
    mr = mr.sort_values('mean_expr', ascending=False).drop_duplicates('Gene_Symbol')
    expr_mr = mr.set_index('Gene_Symbol')[T1_COLS + T2_COLS]
    print(f"mRNA: {len(expr_mr)} unique gene symbols")

    for expr, name, tag in [(expr_mi, 'miRNA', 'DEM'), (expr_mr, 'Gene', 'DEG')]:
        res, info = moderated_paired_analysis(expr, name)
        cand = select_candidates(res)
        wider = select_moderated_only(res)
        res.to_csv(os.path.join(args.out, f'{name}_moderated_all_results.csv'), index=False)
        cand.to_csv(os.path.join(args.out, f'{tag}_candidates.csv'), index=False)
        wider.to_csv(os.path.join(args.out, f'{tag}_moderated_only_candidates.csv'), index=False)

        sens = female_pair_sensitivity(expr)
        overlap = set(cand[name]) & sens
        print(f"\n=== {name} ===")
        print(f"eBayes prior df d0={info['prior_df_d0']:.3f}, s0^2={info['prior_var_s02']:.4f}")
        print(f"{tag} candidates (dual-test P<{P_THR}, |FC|>={2**FC_THR:.1f}, 3-pair consistent): {len(cand)}")
        print(f"  up in T1: {(cand['logFC'] < 0).sum()}, up in T2: {(cand['logFC'] > 0).sum()}")
        print(f"moderated-only candidates (supplementary): {len(wider)}")
        print(f"FDR (BH) < 0.05: {(res['adj.P.Val'] < 0.05).sum()} (min adj.P = {res['adj.P.Val'].min():.4f})")
        print(f"female-pair sensitivity overlap: {len(overlap)}/{len(cand)}"
              f" ({100 * len(overlap) / max(len(cand), 1):.1f}%)")
        print(cand[[name, 'logFC', 't_moderated', 'P.Value', 'adj.P.Val', 'P.ordinary_t']].to_string(index=False))

    print(f"\nResults written to {args.out}")


if __name__ == '__main__':
    main()
