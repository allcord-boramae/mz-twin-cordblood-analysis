"""
Computational cell-type deconvolution of MZ twin cord blood mRNA data.

Input: processed_mRNA_expression_matrix.txt exactly as deposited in
NCBI GEO (GSE326366).

Methods (identical to the manuscript):
 1. NNLS (CIBERSORT-like) regression against a curated binary signature
    matrix of 11 cord blood cell types (LM22-inspired marker genes;
    Newman et al., Nat Methods 2015) -> estimated proportions (%).
 2. Mean marker z-score per cell type (sample-wise z-scored expression).
 3. Simplified ssGSEA score (normalized mean rank of marker genes).
Within-pair (T1 vs T2) comparisons use the exact paired Wilcoxon
signed-rank test for all three methods; note that with n = 3 pairs the
minimum attainable exact two-sided P value is 0.25, so results are
descriptive.

Usage:
  python 02_deconvolution.py --geo-dir <dir> --out ../results
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy import stats

SAMPLES = ['MZ1-T1', 'MZ1-T2', 'MZ2-T1', 'MZ2-T2', 'MZ3-T1', 'MZ3-T2']
PAIRS = ['MZ1', 'MZ2', 'MZ3']

SIGNATURE = {
    'T cells CD4': ['CD3D','CD3E','CD4','IL7R','LEF1','TCF7','CD27','CD28','LCK','ITK'],
    'T cells CD8': ['CD3D','CD3E','CD8A','CD8B','GZMK','EOMES','CCL5','NKG7','PRF1','GZMA'],
    'Treg': ['FOXP3','IL2RA','CTLA4','IKZF2','TNFRSF18','CD3D','CD3E','CD4'],
    'B cells': ['CD19','MS4A1','CD79A','CD79B','BLK','BANK1','PAX5','TCL1A','IGHM','IGHD'],
    'NK cells': ['NCAM1','KLRD1','NKG7','GNLY','GZMB','PRF1','KLRB1','KLRF1','FCGR3A','NCR1'],
    'Monocytes': ['CD14','CD68','LYZ','S100A8','S100A9','FCN1','VCAN','CSF1R','ITGAM','CST3'],
    'Neutrophils': ['FCGR3B','CEACAM8','CSF3R','CXCR1','CXCR2','S100A12','MMP9','CAMP','LCN2','AQP9'],
    'Erythroid': ['GYPA','HBB','HBA1','HBA2','ALAS2','SLC4A1','ANK1','EPB42','GATA1','KLF1'],
    'HSC/Progenitor': ['CD34','KIT','FLT3','HOXA9','MEIS1','GATA2','RUNX1','MPL','THY1','PROM1'],
    'pDC': ['IL3RA','CLEC4C','IRF7','TCF4','LILRA4','JCHAIN','GZMB','ITM2C','PLD4','SERPINF1'],
    'mDC': ['CD1C','ITGAX','FCER1A','CLEC10A','HLA-DQA1','HLA-DPA1','BATF3','IRF4','CD1E','THBD'],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geo-dir', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(os.path.join(args.geo_dir, 'processed_mRNA_expression_matrix.txt'), sep='\t')
    df_gene = df.groupby('Gene_Symbol')[SAMPLES].max().reset_index()

    cell_types = list(SIGNATURE.keys())
    all_sig = set(g for genes in SIGNATURE.values() for g in genes)
    sig_found = sorted(all_sig & set(df_gene['Gene_Symbol']))
    print(f"Signature genes found: {len(sig_found)}/{len(all_sig)}")

    S = np.array([[1.0 if g in SIGNATURE[ct] else 0.0 for ct in cell_types] for g in sig_found])
    expr_sig = df_gene[df_gene['Gene_Symbol'].isin(sig_found)] \
        .set_index('Gene_Symbol').loc[sig_found]

    # Method 1: NNLS proportions
    props = {}
    for sample in SAMPLES:
        b = expr_sig[sample].values
        b_std = (b - b.mean()) / (b.std() + 1e-10)
        S_std = (S - S.mean(axis=0)) / (S.std(axis=0) + 1e-10)
        w, _ = nnls(S_std, b_std)
        props[sample] = dict(zip(cell_types, w / (w.sum() + 1e-10) * 100))
    prop_df = pd.DataFrame(props).T.round(2)
    print("\nNNLS estimated proportions (%):")
    print(prop_df.to_string())

    # Methods 2 & 3: marker z-score and simplified ssGSEA score
    expr_all = df_gene.set_index('Gene_Symbol')[SAMPLES]
    marker, ssgsea = {}, {}
    ranks = {s: stats.rankdata(expr_all[s].values) / len(expr_all) for s in SAMPLES}
    for sample in SAMPLES:
        vals = expr_all[sample].values
        mu, sigma = vals.mean(), vals.std()
        marker[sample], ssgsea[sample] = {}, {}
        for ct, genes in SIGNATURE.items():
            avail = [g for g in genes if g in expr_all.index]
            idx = expr_all.index.get_indexer(avail)
            marker[sample][ct] = float(((expr_all.loc[avail, sample].values - mu) / (sigma + 1e-10)).mean())
            ssgsea[sample][ct] = float(ranks[sample][idx].mean())
    marker_df = pd.DataFrame(marker).T.round(4)
    ssgsea_df = pd.DataFrame(ssgsea).T.round(4)

    def within_pair_stats(score_df, value_label):
        rows = []
        for ct in cell_types:
            t1 = np.array([score_df.loc[f'{p}-T1', ct] for p in PAIRS])
            t2 = np.array([score_df.loc[f'{p}-T2', ct] for p in PAIRS])
            d = t1 - t2
            if np.allclose(d, 0):
                p_w = 1.0
            else:
                p_w = stats.wilcoxon(t1, t2, zero_method='wilcox', mode='exact').pvalue
            rows.append({
                'Cell Type': ct,
                f'Mean T1 ({value_label})': t1.mean(),
                f'Mean T2 ({value_label})': t2.mean(),
                'Mean diff T1-T2': d.mean(),
                'diff_MZ1': d[0], 'diff_MZ2': d[1], 'diff_MZ3': d[2],
                'Wilcoxon exact P': p_w,
                '3-pair consistent': bool(np.all(d > 0) or np.all(d < 0)),
            })
        return pd.DataFrame(rows)

    res_nnls = within_pair_stats(prop_df, '%')
    res_marker = within_pair_stats(marker_df, 'z')
    res_ssgsea = within_pair_stats(ssgsea_df, 'rank')
    print("\nWithin-pair comparison, NNLS proportions (exact paired Wilcoxon signed-rank):")
    print(res_nnls.round(3).to_string(index=False))
    print("\nNote: with n=3 pairs the minimum attainable exact two-sided P is 0.25.")

    prop_df.to_csv(os.path.join(args.out, 'deconvolution_NNLS_proportions.csv'))
    marker_df.to_csv(os.path.join(args.out, 'deconvolution_marker_zscores.csv'))
    ssgsea_df.to_csv(os.path.join(args.out, 'deconvolution_ssgsea_scores.csv'))
    res_nnls.to_csv(os.path.join(args.out, 'deconvolution_within_pair_stats_NNLS.csv'), index=False)
    res_marker.to_csv(os.path.join(args.out, 'deconvolution_within_pair_stats_marker.csv'), index=False)
    res_ssgsea.to_csv(os.path.join(args.out, 'deconvolution_within_pair_stats_ssgsea.csv'), index=False)
    print(f"\nResults written to {args.out}")


if __name__ == '__main__':
    main()
