import matplotlib
matplotlib.use('Agg')
# CRITICAL: SVG fonttype='none' → text as <text> elements, not paths → Illustrator editable
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.font_manager as fm

# ============================================================
# Nature-inspired color palette
# ============================================================
N_BLUE   = '#4472C4'   # steel blue
N_ORANGE = '#ED7D31'   # rust orange
N_GREEN  = '#70AD47'   # meadow green (Nature signature)
N_RED    = '#E74C3C'   # muted red
N_GRAY   = '#BDC3C7'   # silver gray (for unlabeled only)

# ============================================================
# CHART 1: Species Distribution Donut
# ============================================================
species = {
    'Nurse shark\n(Ginglymostoma cirratum)': 226,
    'Bamboo shark\n(Chiloscyllium plagiosum)': 18,
    'Spotted wobbegong\n(Orectolobus maculatus)': 11,
    'Spiny dogfish\n(Squalus acanthias)': 2,
}
colors1 = [N_BLUE, N_ORANGE, N_GREEN, N_RED]
explode1 = (0, 0.06, 0.06, 0.06)

fig, ax = plt.subplots(figsize=(9, 7))

wedges, texts = ax.pie(
    list(species.values()),
    labels=list(species.keys()),
    colors=colors1,
    startangle=90,
    labeldistance=1.15,
    explode=explode1,
    textprops={'fontsize': 11, 'fontfamily': 'Arial', 'linespacing': 1.3},
    wedgeprops={'width': 0.40, 'edgecolor': 'white', 'linewidth': 2.5},
    pctdistance=0.0,
)

total = sum(species.values())
for w, (label, val) in zip(wedges, species.items()):
    ang = (w.theta2 + w.theta1) / 2
    x = 0.55 * np.cos(np.deg2rad(ang))
    y = 0.55 * np.sin(np.deg2rad(ang))
    pct = val / total * 100
    ax.text(x, y, f'{pct:.1f}%\n(n={val})', ha='center', va='center',
            fontsize=10, fontweight='bold', fontfamily='Arial', color='white')

ax.text(0, 0, 'N = 257\n(excl. 4 ENGINEERING)', ha='center', va='center',
        fontsize=12, fontweight='bold', fontfamily='Arial', color='#333333')

ax.set_title('VNAR Sequence Collection: Species Distribution',
             fontsize=15, fontweight='bold', fontfamily='Arial', pad=25, color='#333333')

fig.text(0.5, 0.02,
         'Bamboo shark: 8 SOLUBLE + 10 INSOLUBLE  |  Others: SOLUBLE only  |  NCBI 220 seqs: nurse shark origin',
         ha='center', fontsize=8, fontfamily='Arial', color='#666666', style='italic')

fig.tight_layout(rect=[0, 0.04, 1, 0.96])
out1 = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\figures\Fig_Chart_Species_Donut.svg'
fig.savefig(out1, format='svg', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out1}')

# ============================================================
# CHART 2: Solubility Labels Donut
# ============================================================
labels_data = {'SOLUBLE': 27, 'INSOLUBLE': 10, 'Unlabeled NCBI': 220}
colors2 = [N_GREEN, N_RED, N_GRAY]

fig2, ax2 = plt.subplots(figsize=(8, 6.5))
wedges2, texts2 = ax2.pie(
    list(labels_data.values()),
    labels=None,
    colors=colors2,
    startangle=90,
    wedgeprops={'width': 0.40, 'edgecolor': 'white', 'linewidth': 2.5},
    pctdistance=0.0,
)

total2 = sum(labels_data.values())
for w, (label, val) in zip(wedges2, labels_data.items()):
    ang = (w.theta2 + w.theta1) / 2
    x = 0.55 * np.cos(np.deg2rad(ang))
    y = 0.55 * np.sin(np.deg2rad(ang))
    pct = val / total2 * 100
    ax2.text(x, y, f'{pct:.1f}%\n(n={val})', ha='center', va='center',
            fontsize=11, fontweight='bold', fontfamily='Arial', color='white')

legend_labels = [f'{k} (n={v}, {v/total2*100:.1f}%)' for k, v in labels_data.items()]
ax2.legend(wedges2, legend_labels, loc='lower center', bbox_to_anchor=(0.5, -0.08),
          fontsize=10, ncol=3, frameon=False)

ax2.text(0, 0, 'N = 257', ha='center', va='center',
        fontsize=13, fontweight='bold', fontfamily='Arial', color='#333333')

ax2.set_title('VNAR Sequence Collection: Solubility Labels',
             fontsize=15, fontweight='bold', fontfamily='Arial', pad=25, color='#333333')

fig2.tight_layout(rect=[0, 0.08, 1, 0.96])
out2 = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\figures\Fig_Chart_Labels_Donut.svg'
fig2.savefig(out2, format='svg', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out2}')

# ============================================================
# CHART 3: 3-Set Venn — Nature Blue/Orange/Green
# ============================================================
from matplotlib_venn import venn3

values = (4, 1, 0, 7, 2, 2, 3)

fig3, ax3 = plt.subplots(figsize=(10, 8))
v = venn3(subsets=values,
          set_labels=('Level 1a: All-species\nFisher (9 positions)',
                      'Level 1b: Species-controlled\nFisher (6 positions)',
                      'Level 2: Bamboo-only\nFisher (14 positions)'),
          set_colors=(N_BLUE, N_ORANGE, N_GREEN),
          alpha=0.5,
          ax=ax3)

label_map = {
    '100': 'L1a only (4)\n25, 44, 46, 69',
    '010': 'L1b only (1)\n34',
    '001': 'L2 only (7)\n4,27,42,58,84,\n165,167',
    '110': 'L1a\u2229L1b (none)',
    '101': 'L1a\u2229L2 (2)\n39, 47',
    '011': 'L1b\u2229L2 (2)\n49, 50',
    '111': 'ROBUST (3)\n43(G), 83(I), 94(T)',
}

for rid, text in label_map.items():
    lbl = v.get_label_by_id(rid)
    if lbl:
        lbl.set_text(text)
        lbl.set_fontsize(8)
        lbl.set_fontfamily('Arial')

for sid in ['A', 'B', 'C']:
    lbl = v.get_label_by_id(sid)
    if lbl:
        lbl.set_fontsize(12)
        lbl.set_fontweight('bold')
        lbl.set_fontfamily('Arial')

ax3.set_title('Fisher Significance: Cross-Level Overlap\n(ALL_MSA numbering; Level 2 mapped from Bamboo MSA)',
             fontsize=14, fontweight='bold', fontfamily='Arial', pad=20, color='#333333')

fig3.text(0.5, 0.02,
          'Triple overlap: col43(G), col83(I), col94(T) \u2014 robust across all three analysis levels',
          ha='center', fontsize=9, fontfamily='Arial', color='#C55A11', fontweight='bold',
          bbox=dict(facecolor='#FFF5EB', edgecolor=N_ORANGE, boxstyle='round,pad=0.5'))

fig3.tight_layout(rect=[0, 0.05, 1, 0.96])
out3 = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\figures\Fig_Chart_Venn_Fisher.svg'
fig3.savefig(out3, format='svg', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out3}')

print('\nDone! All SVGs use svg.fonttype=none — text is fully editable in Illustrator.')
