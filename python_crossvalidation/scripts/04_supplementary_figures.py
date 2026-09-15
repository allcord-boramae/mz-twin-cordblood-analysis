"""
S3 Fig: paired within-pair changes in estimated cell-type proportions
(slope / spaghetti plot), one panel per cell type.

Reads the NNLS proportions produced by 02_deconvolution.py and draws, for
each of the 11 cell types, the T1 -> T2 change in each of the three MZ
pairs. This presents the deconvolution comparison descriptively, as
appropriate for n = 3 pairs (minimum attainable exact paired Wilcoxon
two-sided P = 0.25).

Usage: python 04_supplementary_figures.py --results ../results --out ../figures
"""
import argparse
import io
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

MAX_WIDTH_IN = 7.5
PAIRS = ['MZ1', 'MZ2', 'MZ3']
PAIR_COLORS = {'MZ1': '#2166AC', 'MZ2': '#E67E22', 'MZ3': '#27AE60'}


def save_plos_tiff(fig, path_base):
    fig_w = fig.get_size_inches()[0]
    dpi_save = min(300, int(600 * MAX_WIDTH_IN / fig_w))
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi_save, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    img = Image.open(buf)
    dpi_meta = img.width / MAX_WIDTH_IN
    img.save(path_base + '.tif', compression='tiff_lzw', dpi=(dpi_meta, dpi_meta))
    fig.savefig(path_base + '.png', dpi=dpi_save, bbox_inches='tight', facecolor='white')
    print(f"{os.path.basename(path_base)}: {img.width}x{img.height}px, "
          f"{dpi_meta:.0f} dpi at {MAX_WIDTH_IN} in width")


DEG_ALL = ['CISD1', 'DNTTIP2', 'FAM91A1', 'IGKV2-29', 'IGLV1-40', 'LOC100128653',
           'LOC100128751', 'LOC100131174', 'LOC100996579', 'LOC440602', 'MIR339',
           'PLEKHM1P', 'RABGGTB', 'SH2D1A', 'SLC52A3', 'SLED1', 'TOMM40L', 'UBE2Q2L',
           'OR2L3', 'PHC1', 'PRKX', 'R3HDM1', 'SGOL1']
DEG_MISSED_FEMALE = {'OR2L3', 'PHC1', 'PRKX', 'R3HDM1', 'SGOL1'}


def fig_s2(out):
    """S2 Fig: sensitivity analysis (three-pair vs. female-only)."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), gridspec_kw={'width_ratios': [1, 1, 1.6]})

    for ax, counts, labels, ylab, title in [
        (axes[0], [8, 8, 8], ['3-pair\nanalysis\n(n=8)', 'Overlap\n(n=8)', 'Female-only\n(n=8)'],
         'Number of DEMs', '(A) DEM overlap\n(100% concordance)'),
        (axes[1], [23, 18, 18], ['3-pair\nanalysis\n(n=23)', 'Overlap\n(n=18)', 'Female-only\n(n=18)'],
         'Number of DEGs', '(B) DEG overlap\n(78.3% concordance)'),
    ]:
        bars = ax.bar(range(3), counts, color=['#2166AC', '#9E7BB5', '#B2182B'],
                      edgecolor='black', linewidth=0.8, width=0.62)
        for b, v in zip(bars, counts):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.25, str(v), ha='center',
                    fontsize=13, fontweight='bold')
        ax.set_xticks(range(3))
        ax.set_xticklabels(labels, fontsize=10.5, fontweight='bold')
        ax.set_ylabel(ylab, fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
        ax.set_ylim(0, max(counts) * 1.18)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='y', labelsize=10)

    ax = axes[2]
    n = len(DEG_ALL)
    for i, g in enumerate(DEG_ALL):
        y = n - 1 - i
        ax.barh(y, 1, left=0, color='#2166AC', edgecolor='white', height=0.92)
        if g in DEG_MISSED_FEMALE:
            ax.barh(y, 1, left=1.05, color='#F2F2F2', edgecolor='#CCCCCC', height=0.92)
            ax.text(1.575, y, '×', ha='center', va='center', fontsize=11,
                    color='#B2182B', fontweight='bold')
            ax.text(0.5, y, '●', ha='center', va='center', fontsize=8, color='white')
        else:
            ax.barh(y, 1, left=1.05, color='#2166AC', edgecolor='white', height=0.92)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(DEG_ALL)), fontsize=9)
    ax.set_xticks([0.5, 1.575])
    ax.set_xticklabels(['3-pair\n(MZ1+MZ2+MZ3)\ndual-test criteria', 'Female-only\n(MZ1+MZ2)\nFC + direction only'],
                       fontsize=10.5, fontweight='bold')
    ax.set_xlim(-0.02, 2.12)
    ax.set_title('(C) Per-gene concordance', fontsize=13, fontweight='bold', pad=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False)
    import matplotlib.patches as mpatches
    ax.legend(handles=[mpatches.Patch(color='#2166AC', label='Detected'),
                       mpatches.Patch(facecolor='#F2F2F2', edgecolor='#CCCCCC', label='Not detected')],
              loc='upper left', bbox_to_anchor=(0.0, -0.14), ncol=2, fontsize=9.5, frameon=False)
    fig.text(0.5, 0.005,
             'Three-pair analysis: |FC| ≥ 1.5, moderated and ordinary paired t-test P < 0.05, 3-pair '
             'directional consistency. Female-only analysis: |FC| ≥ 1.5 and 2-pair directional concordance '
             '(no significance test is possible with n = 2 pairs).',
             ha='center', fontsize=9, style='italic')
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    save_plos_tiff(fig, os.path.join(out, 'S2_Fig'))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    fig_s2(args.out)

    prop = pd.read_csv(os.path.join(args.results, 'deconvolution_NNLS_proportions.csv'), index_col=0)
    stats = pd.read_csv(os.path.join(args.results, 'deconvolution_within_pair_stats_NNLS.csv'))
    cell_types = list(prop.columns)

    fig, axes = plt.subplots(3, 4, figsize=(13, 10))
    axes = axes.flatten()

    for i, ct in enumerate(cell_types):
        ax = axes[i]
        vals = []
        for p in PAIRS:
            t1 = prop.loc[f'{p}-T1', ct]
            t2 = prop.loc[f'{p}-T2', ct]
            vals.extend([t1, t2])
            ax.plot([0, 1], [t1, t2], marker='o', markersize=7, lw=2.2,
                    color=PAIR_COLORS[p], alpha=0.85, zorder=3)
        row = stats[stats['Cell Type'] == ct].iloc[0]
        consistent = bool(row['3-pair consistent'])
        title = ct + (' *' if consistent else '')
        ax.set_title(title, fontsize=13, fontweight='bold',
                     color='#B2182B' if consistent else '#222222')
        ax.set_xlim(-0.35, 1.35)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['T1', 'T2'], fontsize=12, fontweight='bold')
        vmax = max(vals) if vals else 1
        pad = max(0.6, (vmax - min(vals)) * 0.35) if vals else 1
        ax.set_ylim(max(-0.3, min(vals) - pad), vmax + pad)
        ax.tick_params(axis='y', labelsize=10)
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        d = row['Mean diff T1-T2']
        ax.text(0.5, 0.02, f'mean Δ = {-d:+.2f} pp', transform=ax.transAxes,
                ha='center', fontsize=9.5, style='italic', color='#555555')

    axes[11].axis('off')
    handles = [Line2D([0], [0], color=PAIR_COLORS[p], marker='o', lw=2.2, label=p)
               for p in PAIRS]
    handles.append(Line2D([0], [0], color='none', label='* 3-pair consistent direction'))
    axes[11].legend(handles=handles, loc='center', fontsize=12, frameon=True,
                    title='Twin pair', title_fontsize=13)

    fig.text(0.06, 0.5, 'Estimated cell-type proportion (%)', va='center',
             rotation='vertical', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0.07, 0, 1, 1])
    save_plos_tiff(fig, os.path.join(args.out, 'S3_Fig'))
    plt.close(fig)


if __name__ == '__main__':
    main()
