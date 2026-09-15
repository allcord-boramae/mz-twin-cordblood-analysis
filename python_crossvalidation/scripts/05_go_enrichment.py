"""
GO enrichment analysis of the 23 candidate DEGs.

Inputs:
  - processed_mRNA_expression_matrix.txt (GEO GSE326366) -> background universe
    (all unique gene symbols on the Human Gene 2.0 ST array, n = 29,383)
  - annotation/HuGene20ST_GO_annotation.tsv.gz (gene -> GO term mapping from
    the Affymetrix Human Gene 2.0 ST array annotation; BP/CC/MF)
  - results/DEG_candidates.csv (the 23 candidate DEGs)

Method: one-sided Fisher's exact test (overrepresentation) per GO term.
As in standard enrichment tools (e.g., clusterProfiler, DAVID), only terms
annotated to at least one candidate DEG are testable; Benjamini-Hochberg
FDR correction is applied across all such tested terms (the number of
tested terms is reported in the output and in the manuscript Methods).
Full per-term results (raw P, adjusted P, term size, hit genes) are
written to results/GO_enrichment_full_results.csv.

Usage: python 05_go_enrichment.py --geo-dir <dir> --annotation
       ../annotation/HuGene20ST_GO_annotation.tsv.gz
       --deg ../results/DEG_candidates.csv --out ../results
"""
import argparse
import os
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--geo-dir', required=True)
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--deg', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    mr = pd.read_csv(os.path.join(args.geo_dir, 'processed_mRNA_expression_matrix.txt'), sep='\t')
    universe = set(str(g).strip() for g in mr['Gene_Symbol'].dropna().unique())
    universe.discard('---')
    N = len(universe)
    print(f"Background universe: {N} unique gene symbols")

    ann = pd.read_csv(args.annotation, sep='\t')
    go2genes = defaultdict(set)
    go_meta = {}
    for g, gid, term, ont in ann.itertuples(index=False):
        if g in universe:
            go2genes[gid].add(g)
            go_meta[gid] = (term, ont)
    print(f"GO terms with >=1 annotated gene in universe: {len(go2genes)}")

    deg = set(pd.read_csv(args.deg)['Gene'])
    n_deg = len(deg)
    print(f"Candidate DEGs: {n_deg}")

    rows = []
    for gid, genes in go2genes.items():
        overlap = deg & genes
        k = len(overlap)
        if k == 0:
            continue
        K = len(genes)
        table = [[k, n_deg - k], [K - k, N - K - n_deg + k]]
        odds, p = fisher_exact(table, alternative='greater')
        term, ont = go_meta[gid]
        rows.append({'Ontology': ont, 'GO_ID': gid, 'GO_Term': term,
                     'K_universe': K, 'k_DEG': k, 'P_value': p,
                     'odds_ratio': odds, 'genes': ','.join(sorted(overlap))})

    res = pd.DataFrame(rows).sort_values('P_value').reset_index(drop=True)
    m = len(res)
    print(f"Tested GO terms (>=1 candidate DEG annotated): m = {m}")

    ranks = np.arange(1, m + 1)
    adj = res['P_value'].to_numpy() * m / ranks
    for i in range(m - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    res['FDR_BH'] = np.minimum(adj, 1.0)

    sig = res[res['FDR_BH'] < 0.05]
    print(f"FDR < 0.05: {len(sig)} terms "
          f"(BP {sum(sig.Ontology == 'BP')}, CC {sum(sig.Ontology == 'CC')}, MF {sum(sig.Ontology == 'MF')})")
    print(sig[['Ontology', 'GO_ID', 'GO_Term', 'K_universe', 'k_DEG', 'P_value', 'FDR_BH']]
          .head(30).to_string(index=False))

    res.to_csv(os.path.join(args.out, 'GO_enrichment_full_results.csv'), index=False)
    print(f"\nFull results ({m} tested terms) written to GO_enrichment_full_results.csv")


if __name__ == '__main__':
    main()
