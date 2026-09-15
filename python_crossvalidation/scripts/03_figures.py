"""
Generate main figures (Fig 1-4) for the manuscript.

Revision R1 changes (reviewer #1, minor comments 2-3):
 - all figures exported as LZW-compressed TIFF, physical width 7.5 in,
   effective resolution 300-600 dpi (PLOS ONE figure specifications)
 - Figures 2, 3 and 4: all text and number font sizes increased by >= 50%
   relative to the originally submitted versions

Figures are drawn from the candidate DEM/DEG statistics produced by
01_differential_expression.py (values embedded below for self-containment;
they are identical to results/DEM_candidates.csv and DEG_candidates.csv).

Usage:  python 03_figures.py --out ../figures
"""
import argparse
import io
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from PIL import Image

MAX_WIDTH_IN = 7.5     # PLOS ONE maximum figure width
FS2 = 1.5              # font scale for Figures 2-4 (reviewer: >= +50%)

# ---------------------------------------------------------------------------
# Data (from results/DEM_candidates.csv / DEG_candidates.csv)
# ---------------------------------------------------------------------------
DEM_DATA = [
    ("hsa-miR-6779-5p",  -0.69, -0.64, -0.71),
    ("hsa-miR-1292-5p",  -1.62, -1.58, -1.79),
    ("hsa-miR-6891-5p",  -0.81, -0.89, -0.69),
    ("hsa-miR-4253",     -0.79, -0.91, -0.71),
    ("hsa-miR-128-3p",   -1.47, -0.84, -1.21),
    ("hsa-miR-148b-3p",  -0.62, -1.01, -1.10),
    ("hsa-miR-3679-5p",  -0.87, -0.54, -0.49),
    ("hsa-miR-328-3p",   -0.62, -1.48, -1.18),
]
DEG_DATA = [
    ("LOC100131174", +0.78, +0.77, +0.70, "T2"),
    ("CISD1",        +0.65, +0.71, +0.48, "T2"),
    ("UBE2Q2L",      +0.67, +0.95, +0.66, "T2"),
    ("LOC440602",    +0.77, +0.49, +0.60, "T2"),
    ("LOC100996579", +0.76, +0.41, +0.59, "T2"),
    ("IGLV1-40",     +1.33, +0.69, +1.07, "T2"),
    ("SLC52A3",      +0.60, +0.87, +0.46, "T2"),
    ("MIR339",       +0.50, +0.96, +0.59, "T2"),
    ("LOC100128653", +0.63, +0.91, +0.39, "T2"),
    ("SLED1",        -0.79, -0.64, -0.71, "T1"),
    ("SH2D1A",       -0.62, -0.56, -0.73, "T1"),
    ("IGKV2-29",     -0.75, -0.57, -0.60, "T1"),
    ("FAM91A1",      -0.98, -0.78, -0.69, "T1"),
    ("TOMM40L",      -0.66, -0.59, -0.86, "T1"),
    ("LOC100128751", -1.13, -0.77, -0.82, "T1"),
    ("PRKX",         -0.70, -0.44, -0.62, "T1"),
    ("DNTTIP2",      -0.55, -0.82, -0.53, "T1"),
    ("PHC1",         -0.50, -0.66, -0.89, "T1"),
    ("R3HDM1",       -0.43, -0.60, -0.81, "T1"),
    ("PLEKHM1P",     -0.72, -0.87, -1.38, "T1"),
    ("SGOL1",        -0.70, -0.37, -0.80, "T1"),
    ("OR2L3",        -0.36, -0.72, -0.75, "T1"),
    ("RABGGTB",      -0.84, -0.83, -0.39, "T1"),
]
# 2026-09-13 실DB 재검증 반영: 17쌍 (328->TOMM40L 삭제: miRDB v6 부재,
# 4253->PRKX는 TS 단독으로 강등: miRTarBase v9에 CLIP-seq 기록 부재).
# 근거: miRDB_v6.0_prediction_result.txt.gz, miRTarBase_MTI.xlsx,
# TargetScan8.0 per-miRNA exports (클로드챗 폴더, 2026-03-22).
INTERACTIONS = [
    ("hsa-miR-128-3p",  "CISD1",   "Anti-corr", ["miRTarBase"], "CLASH"),
    ("hsa-miR-1292-5p", "CISD1",   "Anti-corr", ["TargetScan"], ""),
    ("hsa-miR-6779-5p", "CISD1",   "Anti-corr", ["TargetScan"], ""),
    ("hsa-miR-6891-5p", "CISD1",   "Anti-corr", ["TargetScan"], ""),
    ("hsa-miR-128-3p",  "PRKX",    "Same", ["miRDB", "TargetScan"], ""),
    ("hsa-miR-4253",    "PRKX",    "Same", ["TargetScan"], ""),
    ("hsa-miR-6779-5p", "PHC1",    "Same", ["TargetScan", "miRTarBase"], "PAR-CLIP"),
    ("hsa-miR-1292-5p", "PHC1",    "Same", ["TargetScan"], ""),
    ("hsa-miR-128-3p",  "TOMM40L", "Same", ["miRTarBase"], "PAR-CLIP"),
    ("hsa-miR-6779-5p", "TOMM40L", "Same", ["TargetScan"], ""),
    ("hsa-miR-6891-5p", "TOMM40L", "Same", ["TargetScan"], ""),
    ("hsa-miR-4253",    "TOMM40L", "Same", ["TargetScan"], ""),
    ("hsa-miR-6779-5p", "SGOL1",   "Same", ["TargetScan"], ""),
    ("hsa-miR-6891-5p", "SGOL1",   "Same", ["TargetScan"], ""),
    ("hsa-miR-6891-5p", "SH2D1A",  "Same", ["TargetScan"], ""),
    ("hsa-miR-6779-5p", "DNTTIP2", "Same", ["TargetScan"], ""),
    ("hsa-miR-4253",    "DNTTIP2", "Same", ["TargetScan"], ""),
]
# (Ontology, GO ID, term, K universe, k DEG, P, FDR, genes)
GO_TERMS = [
    ('BP', 'GO:2000696', 'regulation of epithelial cell differentiation involved in kidney development', 1, 1, 0.0007827655447027193, 0.04372015863287523, 'PRKX'),
    ('BP', 'GO:0060562', 'epithelial tube morphogenesis', 1, 1, 0.0007827655447027193, 0.04372015863287523, 'PRKX'),
    ('CC', 'GO:0030892', 'mitotic cohesin complex', 2, 1, 0.001564944987643017, 0.04372015863287523, 'SGOL1'),
    ('BP', 'GO:0032218', 'riboflavin transport', 3, 1, 0.002346538747735738, 0.04372015863287523, 'SLC52A3'),
    ('MF', 'GO:0032217', 'riboflavin transporter activity', 3, 1, 0.002346538747735738, 0.04372015863287523, 'SLC52A3'),
    ('CC', 'GO:0005968', 'Rab-protein geranylgeranyltransferase complex', 4, 1, 0.003127547243610553, 0.04372015863287523, 'RABGGTB'),
    ('MF', 'GO:0004663', 'Rab geranylgeranyltransferase activity', 4, 1, 0.003127547243610553, 0.04372015863287523, 'RABGGTB'),
    ('BP', 'GO:0045132', 'meiotic chromosome segregation', 4, 1, 0.003127547243610553, 0.04372015863287523, 'SGOL1'),
    ('BP', 'GO:0006771', 'riboflavin metabolic process', 5, 1, 0.003907970893612152, 0.04372015863287523, 'SLC52A3'),
    ('BP', 'GO:0018344', 'protein geranylgeranylation', 6, 1, 0.004687810115800419, 0.04372015863287523, 'RABGGTB'),
    ('MF', 'GO:0015288', 'porin activity', 6, 1, 0.004687810115800419, 0.04372015863287523, 'TOMM40L'),
    ('BP', 'GO:0010457', 'centriole-centriole cohesion', 6, 1, 0.004687810115800419, 0.04372015863287523, 'SGOL1'),
    ('MF', 'GO:0004691', 'cAMP-dependent protein kinase activity', 7, 1, 0.005467065327950622, 0.04372015863287523, 'PRKX'),
    ('CC', 'GO:0001739', 'sex chromatin', 7, 1, 0.005467065327950622, 0.04372015863287523, 'PHC1'),
    ('CC', 'GO:0000779', 'condensed chromosome, centromeric region', 7, 1, 0.005467065327950622, 0.04372015863287523, 'SGOL1'),
    ('CC', 'GO:0005741', 'mitochondrial outer membrane', 151, 2, 0.006183705643952328, 0.04372015863287523, 'CISD1,TOMM40L'),
    ('BP', 'GO:0043457', 'regulation of cellular respiration', 8, 1, 0.006245736947553603, 0.04372015863287523, 'CISD1'),
    ('BP', 'GO:0008608', 'attachment of spindle microtubules to kinetochore', 8, 1, 0.006245736947553603, 0.04372015863287523, 'SGOL1'),
    ('CC', 'GO:0046930', 'pore complex', 9, 1, 0.007023825391815946, 0.04657905259835838, 'TOMM40L'),
    ('CC', 'GO:0005742', 'mitochondrial outer membrane translocase complex', 10, 1, 0.007801331077660165, 0.04680798646596099, 'TOMM40L'),
    ('MF', 'GO:0004659', 'prenyltransferase activity', 10, 1, 0.007801331077660165, 0.04680798646596099, 'RABGGTB'),
    ('BP', 'GO:0060993', 'kidney morphogenesis', 11, 1, 0.008578254421724918, 0.04699391552771042, 'PRKX'),
    ('BP', 'GO:0016574', 'histone ubiquitination', 11, 1, 0.008578254421724918, 0.04699391552771042, 'PHC1'),
    ('CC', 'GO:0035102', 'PRC1 complex', 12, 1, 0.009354595840365128, 0.04714716303544025, 'PHC1'),
    ('MF', 'GO:0015266', 'protein channel activity', 12, 1, 0.009354595840365128, 0.04714716303544025, 'TOMM40L'),
    ('CC', 'GO:0000780', 'condensed nuclear chromosome, centromeric region', 13, 1, 0.01013035574965224, 0.04909326247908392, 'SGOL1'),
]


def save_plos_tiff(fig, path_base):
    """Save as LZW TIFF at physical width 7.5 in with 300-600 dpi, plus PNG."""
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


# ---------------------------------------------------------------------------
# Figure 1: heatmap of within-pair differences (resolution improved only)
# ---------------------------------------------------------------------------
def fig1(out):
    fig, (ax_dem, ax_deg) = plt.subplots(1, 2, figsize=(10, 9),
        gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.5})

    dem_names = [d[0] for d in DEM_DATA]
    dem_vals = np.array([[d[1], d[2], d[3]] for d in DEM_DATA])
    im1 = ax_dem.imshow(dem_vals, cmap='RdBu_r', aspect='auto', vmin=-2.5, vmax=2.5)
    ax_dem.set_xticks([0, 1, 2])
    ax_dem.set_xticklabels(['MZ1', 'MZ2', 'MZ3'], fontsize=11, fontweight='bold')
    ax_dem.set_yticks(range(len(dem_names)))
    ax_dem.set_yticklabels(dem_names, fontsize=10, style='italic')
    ax_dem.set_title('A. DEMs (n = 8)\nAll T1-upregulated', fontsize=13, fontweight='bold', pad=12)
    for i in range(len(dem_names)):
        for j in range(3):
            color = 'white' if abs(dem_vals[i, j]) > 1.2 else 'black'
            ax_dem.text(j, i, f'{dem_vals[i, j]:.2f}', ha='center', va='center',
                        fontsize=9, color=color, fontweight='bold')

    deg_names = [d[0] for d in DEG_DATA]
    deg_vals = np.array([[d[1], d[2], d[3]] for d in DEG_DATA])
    im2 = ax_deg.imshow(deg_vals, cmap='RdBu_r', aspect='auto', vmin=-1.5, vmax=1.5)
    ax_deg.set_xticks([0, 1, 2])
    ax_deg.set_xticklabels(['MZ1', 'MZ2', 'MZ3'], fontsize=11, fontweight='bold')
    ax_deg.set_yticks(range(len(deg_names)))
    ax_deg.set_yticklabels(deg_names, fontsize=9)
    ax_deg.set_title('B. DEGs (n = 23)\n14 T1↑ / 9 T2↑', fontsize=13, fontweight='bold', pad=12)
    key_genes = {'CISD1', 'TOMM40L', 'PRKX', 'PHC1', 'SGOL1'}
    for i, name in enumerate(deg_names):
        if name in key_genes:
            ax_deg.get_yticklabels()[i].set_fontweight('bold')
            ax_deg.get_yticklabels()[i].set_color('#C00000')
    for i in range(len(deg_names)):
        for j in range(3):
            color = 'white' if abs(deg_vals[i, j]) > 0.9 else 'black'
            ax_deg.text(j, i, f'{deg_vals[i, j]:+.2f}', ha='center', va='center',
                        fontsize=8, color=color)
    ax_deg.axhline(y=8.5, color='black', linewidth=2)
    ax_deg.text(2.7, 4, 'T2↑', fontsize=11, fontweight='bold', color='#2166AC', va='center')
    ax_deg.text(2.7, 16, 'T1↑', fontsize=11, fontweight='bold', color='#B2182B', va='center')

    cb1 = fig.colorbar(im1, ax=ax_dem, shrink=0.4, pad=0.02, aspect=15)
    cb1.set_label('log₂(T2/T1)', fontsize=10)
    cb2 = fig.colorbar(im2, ax=ax_deg, shrink=0.4, pad=0.02, aspect=15)
    cb2.set_label('log₂(T2/T1)', fontsize=10)
    plt.tight_layout()
    save_plos_tiff(fig, os.path.join(out, 'Fig1'))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2: miRNA-mRNA network (fonts x1.5)
# ---------------------------------------------------------------------------
def fig2(out):
    f = FS2
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(-1, 11)
    ax.set_ylim(-1.3, 10)
    ax.axis('off')

    mirnas = sorted(set(i[0] for i in INTERACTIONS))
    genes = sorted(set(i[1] for i in INTERACTIONS))
    mir_x, gene_x = 1.5, 8.5
    mir_y = {m: 9 - i * (8.5 / (len(mirnas) - 1)) for i, m in enumerate(mirnas)}
    gene_y = {g: 9 - i * (8.5 / (len(genes) - 1)) for i, g in enumerate(genes)}

    for m, y in mir_y.items():
        ax.add_patch(FancyBboxPatch((mir_x - 1.5, y - 0.26), 3.0, 0.52,
                     boxstyle="round,pad=0.1", facecolor='#2166AC',
                     edgecolor='#14396A', linewidth=1.5, alpha=0.9))
        ax.text(mir_x, y, m.replace('hsa-', ''), ha='center', va='center',
                fontsize=8 * f, color='white', fontweight='bold')

    gene_colors = {"CISD1": '#E74C3C', "TOMM40L": '#E67E22', "PRKX": '#27AE60',
                   "PHC1": '#8E44AD', "SGOL1": '#2980B9'}
    for g, y in gene_y.items():
        color = gene_colors.get(g, '#95A5A6')
        ax.add_patch(FancyBboxPatch((gene_x - 1.1, y - 0.26), 2.2, 0.52,
                     boxstyle="round,pad=0.1", facecolor=color,
                     edgecolor='#333333', linewidth=1.5, alpha=0.9))
        ax.text(gene_x, y, g, ha='center', va='center',
                fontsize=9 * f, color='white', fontweight='bold')

    for mirna, gene, direction, dbs, _ in INTERACTIONS:
        y1, y2 = mir_y[mirna], gene_y[gene]
        n_db = len(dbs)
        if direction == "Anti-corr":
            color, style, lw = '#E74C3C', '-', (2.5 if n_db >= 2 else 1.8)
        elif n_db >= 2:
            color, style, lw = '#27AE60', '-', 2.5
        else:
            color, style, lw = '#BDC3C7', '--', 1.0
        ax.annotate('', xy=(gene_x - 1.15, y2), xytext=(mir_x + 1.55, y1),
                    arrowprops=dict(arrowstyle='->', color=color, linestyle=style,
                                    lw=lw, connectionstyle='arc3,rad=0.05'))
        if n_db >= 2:
            abbr = {'miRDB': 'miRDB', 'TargetScan': 'TS', 'miRTarBase': 'MTB'}
            mid_x = (mir_x + 1.55 + gene_x - 1.15) / 2
            mid_y = (y1 + y2) / 2
            ax.text(mid_x, mid_y + 0.15, '+'.join(abbr[d] for d in dbs),
                    fontsize=6 * f, ha='center', color=color, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                              edgecolor=color, alpha=0.8))

    ax.text(mir_x, 9.8, 'DEMs (all T1↑)', ha='center', fontsize=12 * f,
            fontweight='bold', color='#2166AC')
    ax.text(gene_x, 9.8, 'DEG Targets', ha='center', fontsize=12 * f,
            fontweight='bold', color='#333333')

    legend_items = [
        mpatches.Patch(color='#E74C3C', label='Anti-correlated (→ CISD1)'),
        mpatches.Patch(color='#27AE60', label='2-DB validated'),
        mpatches.Patch(color='#BDC3C7', label='1-DB supported'),
        mpatches.Patch(color='#E67E22', label='TOMM40L (4 DEM targets)'),
        mpatches.Patch(color='#8E44AD', label='PHC1 (epigenetic)'),
    ]
    ax.legend(handles=legend_items, loc='lower center', ncol=3, fontsize=8 * f,
              frameon=True, fancybox=True, bbox_to_anchor=(0.5, -0.03))
    plt.tight_layout()
    save_plos_tiff(fig, os.path.join(out, 'Fig2'))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 3: GO enrichment dot plot (fonts x1.5)
# ---------------------------------------------------------------------------
def fig3(out):
    f = FS2
    terms = sorted(GO_TERMS, key=lambda t: t[5])   # by P ascending
    terms = terms[::-1]                            # bottom = most significant
    labels, neglogp, counts, colors = [], [], [], []
    color_map = {'BP': '#E74C3C', 'CC': '#3498DB', 'MF': '#2ECC71'}
    for ont, _, term, K, k, p, fdr, genes in terms:
        name = term if len(term) <= 45 else term[:42] + '...'
        labels.append(name)
        neglogp.append(-np.log10(p))
        counts.append(k)
        colors.append(color_map[ont])

    fig, ax = plt.subplots(figsize=(12, 11))
    ax.scatter(neglogp, range(len(labels)), s=[c * 220 for c in counts],
               c=colors, alpha=0.8, edgecolors='black', linewidth=0.5, zorder=3)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9 * f)
    ax.tick_params(axis='x', labelsize=9 * f)
    ax.set_xlabel('−log₁₀(P-value)', fontsize=12 * f, fontweight='bold')
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    cat = [mpatches.Patch(color='#E74C3C', label='Biological Process'),
           mpatches.Patch(color='#3498DB', label='Cellular Component'),
           mpatches.Patch(color='#2ECC71', label='Molecular Function')]
    legend1 = ax.legend(handles=cat, loc='lower right', fontsize=9 * f,
                        title='GO Category', title_fontsize=10 * f, framealpha=0.9)
    ax.add_artist(legend1)
    size_items = [Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                         markersize=np.sqrt(n * 220 / np.pi) * 0.9, label=f'{n} gene{"s" if n > 1 else ""}',
                         markeredgecolor='black', markeredgewidth=0.5) for n in (1, 2)]
    ax.legend(handles=size_items, loc='center right', fontsize=9 * f,
              title='Gene Count', title_fontsize=10 * f, framealpha=0.9,
              bbox_to_anchor=(1.0, 0.32))
    plt.tight_layout()
    save_plos_tiff(fig, os.path.join(out, 'Fig3'))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: integrated DEM -> DEG -> GO map (fonts x1.5)
# ---------------------------------------------------------------------------
def fig4(out):
    f = FS2
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_xlim(0, 16)
    ax.set_ylim(-0.9, 10.4)
    ax.axis('off')

    col1_x, col2_x, col3_x = 2, 7.5, 13
    key_dems = ["miR-128-3p", "miR-6779-5p", "miR-4253", "miR-1292-5p",
                "miR-6891-5p", "miR-328-3p", "miR-148b-3p", "miR-3679-5p"]
    dem_y = {d: 9 - i * 1.1 for i, d in enumerate(key_dems)}
    deg_y = {"CISD1": 8, "TOMM40L": 6.2, "PRKX": 4.4, "PHC1": 2.6, "SGOL1": 0.8}
    go_terms = [
        ("Mitochondrial\nouter membrane", 8.3, '#E74C3C'),
        ("Regulation of\ncellular respiration", 7.3, '#E74C3C'),
        ("Pore complex /\nTOM complex", 6.3, '#E67E22'),
        ("Porin activity", 5.5, '#E67E22'),
        ("Kidney\ndevelopment", 4.4, '#27AE60'),
        ("cAMP-dependent\nkinase activity", 3.4, '#27AE60'),
        ("PRC1 complex", 2.6, '#8E44AD'),
        ("Histone\nubiquitination", 1.7, '#8E44AD'),
        ("Chromosome\nsegregation", 0.8, '#2980B9'),
    ]

    for d, y in dem_y.items():
        ax.add_patch(FancyBboxPatch((col1_x - 1.35, y - 0.24), 2.7, 0.48,
                     boxstyle="round,pad=0.08", facecolor='#2166AC',
                     edgecolor='#14396A', linewidth=1.2, alpha=0.85))
        ax.text(col1_x, y, d, ha='center', va='center',
                fontsize=7.5 * f, color='white', fontweight='bold')

    deg_colors = {"CISD1": '#E74C3C', "TOMM40L": '#E67E22', "PRKX": '#27AE60',
                  "PHC1": '#8E44AD', "SGOL1": '#2980B9'}
    for g, y in deg_y.items():
        ax.add_patch(FancyBboxPatch((col2_x - 1.0, y - 0.3), 2.0, 0.6,
                     boxstyle="round,pad=0.1", facecolor=deg_colors[g],
                     edgecolor='#333333', linewidth=1.5, alpha=0.9))
        ax.text(col2_x, y, g, ha='center', va='center',
                fontsize=10 * f, color='white', fontweight='bold')

    for term, y, color in go_terms:
        ax.add_patch(FancyBboxPatch((col3_x - 1.6, y - 0.36), 3.2, 0.72,
                     boxstyle="round,pad=0.1", facecolor=color,
                     edgecolor='#333333', linewidth=1, alpha=0.25))
        ax.text(col3_x, y, term, ha='center', va='center',
                fontsize=7.5 * f, color='#222222', fontweight='bold')

    dem_deg_links = [
        ("miR-128-3p", "CISD1", '#E74C3C', 2.0, '-'),
        ("miR-1292-5p", "CISD1", '#E74C3C', 1.5, '-'),
        ("miR-6779-5p", "CISD1", '#E74C3C', 1.5, '-'),
        ("miR-6891-5p", "CISD1", '#E74C3C', 1.5, '-'),
        ("miR-128-3p", "TOMM40L", '#E67E22', 1.2, '--'),
        ("miR-6779-5p", "TOMM40L", '#E67E22', 1.2, '--'),
        ("miR-6891-5p", "TOMM40L", '#E67E22', 1.2, '--'),
        ("miR-4253", "TOMM40L", '#E67E22', 1.2, '--'),
        ("miR-128-3p", "PRKX", '#27AE60', 2.0, '-'),
        ("miR-4253", "PRKX", '#27AE60', 2.0, '-'),
        ("miR-6779-5p", "PHC1", '#8E44AD', 2.0, '-'),
        ("miR-1292-5p", "PHC1", '#8E44AD', 1.2, '--'),
        ("miR-6779-5p", "SGOL1", '#2980B9', 1.2, '--'),
        ("miR-6891-5p", "SGOL1", '#2980B9', 1.2, '--'),
    ]
    for dem, deg, color, lw, style in dem_deg_links:
        ax.annotate('', xy=(col2_x - 1.05, deg_y[deg]), xytext=(col1_x + 1.4, dem_y[dem]),
                    arrowprops=dict(arrowstyle='->', color=color, linestyle=style,
                                    lw=lw, alpha=0.7, connectionstyle='arc3,rad=0.03'))
    deg_go_links = [
        ("CISD1", 8.3, '#E74C3C'), ("CISD1", 7.3, '#E74C3C'),
        ("TOMM40L", 8.3, '#E67E22'), ("TOMM40L", 6.3, '#E67E22'), ("TOMM40L", 5.5, '#E67E22'),
        ("PRKX", 4.4, '#27AE60'), ("PRKX", 3.4, '#27AE60'),
        ("PHC1", 2.6, '#8E44AD'), ("PHC1", 1.7, '#8E44AD'),
        ("SGOL1", 0.8, '#2980B9'),
    ]
    for deg, go_yv, color in deg_go_links:
        ax.annotate('', xy=(col3_x - 1.65, go_yv), xytext=(col2_x + 1.05, deg_y[deg]),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.5, alpha=0.5,
                                    connectionstyle='arc3,rad=0.02'))

    ax.text(col1_x, 10.1, 'DEMs (n = 8)', ha='center', fontsize=13 * f,
            fontweight='bold', color='#2166AC')
    ax.text(col1_x, 9.55, 'All T1-upregulated', ha='center', fontsize=9 * f,
            color='#2166AC', style='italic')
    ax.text(col2_x, 10.1, 'Key DEGs', ha='center', fontsize=13 * f,
            fontweight='bold', color='#333333')
    ax.text(col2_x, 9.55, 'Cross-validated targets', ha='center', fontsize=9 * f,
            color='#555555', style='italic')
    ax.text(col3_x, 10.1, 'GO Enrichment', ha='center', fontsize=13 * f,
            fontweight='bold', color='#333333')
    ax.text(col3_x, 9.55, 'FDR < 0.05', ha='center', fontsize=9 * f,
            color='#555555', style='italic')

    ax.annotate('4 anti-correlated\npairs converge', xy=(col2_x - 1.0, 8),
                xytext=(col2_x - 3.2, 9.3), fontsize=7 * f, color='#E74C3C',
                fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#E74C3C', lw=1),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FDECEA', edgecolor='#E74C3C'))
    ax.annotate('4/8 DEMs\ntarget this gene', xy=(col2_x - 1.0, 6.2),
                xytext=(col2_x - 3.2, 4.9), fontsize=7 * f, color='#E67E22',
                fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#E67E22', lw=1),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FEF5E7', edgecolor='#E67E22'))

    legend_items = [
        Line2D([0], [0], color='#E74C3C', lw=2, label='Anti-correlated (miRNA↑ → mRNA↓)'),
        Line2D([0], [0], color='#27AE60', lw=2, label='2-DB validated (solid)'),
        Line2D([0], [0], color='#BDC3C7', lw=1.2, ls='--', label='1-DB supported (dashed)'),
    ]
    ax.legend(handles=legend_items, loc='lower center', ncol=3, fontsize=8 * f,
              frameon=True, fancybox=True, bbox_to_anchor=(0.5, -0.06))
    plt.tight_layout()
    save_plos_tiff(fig, os.path.join(out, 'Fig4'))
    plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    fig1(args.out)
    fig2(args.out)
    fig3(args.out)
    fig4(args.out)
    print("All figures generated.")

