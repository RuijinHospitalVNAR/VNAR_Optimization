import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# Chart 1: Species Distribution Donut
# ============================================================
species = {
    'Nurse shark\n(6 annotated + 220 unannotated)': 226,
    'Bamboo shark\n(8 soluble + 10 insoluble)': 18,
    'Spotted wobbegong\n(11 soluble)': 11,
    'Spiny dogfish\n(2 soluble)': 2,
}
colors = ['#2F5496', '#F4A261', '#3B7A9E', '#E76F51']
explode = (0, 0.08, 0.08, 0.08)

fig, ax = plt.subplots(figsize=(8, 6))
wedges, texts, autotexts = ax.pie(
    species.values(),
    labels=species.keys(),
    colors=colors,
    autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*sum(species.values()))})',
    startangle=90,
    pctdistance=0.78,
    labeldistance=1.08,
    explode=explode,
    textprops={'fontsize': 11, 'fontfamily': 'Arial'},
    wedgeprops={'width': 0.45, 'edgecolor': 'white', 'linewidth': 2},
)

for at in autotexts:
    at.set_fontsize(10)
    at.set_fontweight('bold')
    at.set_fontfamily('Arial')

ax.set_title('VNAR Sequence Collection: Species Distribution\n(N = 257, 4 ENGINEERING sequences excluded)',
             fontsize=14, fontweight='bold', fontfamily='Arial', pad=25, color='#2F5496')

fig.tight_layout()
out = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\figures\Fig_Chart_Species_Donut.svg'
fig.savefig(out, format='svg', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out}')

# ============================================================
# Chart 2: Solubility Labels Donut
# ============================================================
labels_data = {
    'SOLUBLE (27)': 27,
    'INSOLUBLE (10)': 10,
    'Unlabeled NCBI (220)': 220,
}
colors2 = ['#27AE60', '#E74C3C', '#95A5A6']

fig2, ax2 = plt.subplots(figsize=(8, 6))
wedges2, texts2, autotexts2 = ax2.pie(
    labels_data.values(),
    labels=labels_data.keys(),
    colors=colors2,
    autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*sum(labels_data.values()))})',
    startangle=90,
    pctdistance=0.78,
    labeldistance=1.08,
    textprops={'fontsize': 12, 'fontfamily': 'Arial'},
    wedgeprops={'width': 0.45, 'edgecolor': 'white', 'linewidth': 2},
)

for at in autotexts2:
    at.set_fontsize(11)
    at.set_fontweight('bold')
    at.set_fontfamily('Arial')

ax2.set_title('VNAR Sequence Collection: Solubility Labels\n(N = 257, 4 ENGINEERING sequences excluded)',
              fontsize=14, fontweight='bold', fontfamily='Arial', pad=25, color='#2F5496')

fig2.tight_layout()
out2 = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\figures\Fig_Chart_Labels_Donut.svg'
fig2.savefig(out2, format='svg', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out2}')

# ============================================================
# Chart 3: Venn diagram (Fisher overlap)
# ============================================================
from matplotlib_venn import venn2

fig3, ax3 = plt.subplots(figsize=(7, 6))
v = venn2(
    subsets=(6, 3, 3),
    set_labels=('All-species Fisher\n(9 significant)',
                'Species-controlled\n(6 significant)'),
    set_colors=('#2F5496', '#F4A261'),
    alpha=0.7,
    ax=ax3,
)

# Label the regions
v.get_label_by_id('10').set_text('6 positions\n25, 39, 44, 46, 47, 69')
v.get_label_by_id('01').set_text('3 positions\n34, 49, 50')
v.get_label_by_id('11').set_text('ROBUST (3)\n43(G)  83(I)  94(T)')

for lbl in v.labels:
    if lbl: lbl.set_fontfamily('Arial')

for lbl in v.subset_labels:
    if lbl:
        lbl.set_fontsize(9)
        lbl.set_fontfamily('Arial')

for lbl in v.set_labels:
    if lbl:
        lbl.set_fontsize(12)
        lbl.set_fontweight('bold')
        lbl.set_fontfamily('Arial')

ax3.set_title('Fisher Significance: Cross-Level Overlap\n(ALL_MSA numbering)',
              fontsize=14, fontweight='bold', fontfamily='Arial', pad=20, color='#2F5496')

fig3.tight_layout()
out3 = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\figures\Fig_Chart_Venn_Fisher.svg'
fig3.savefig(out3, format='svg', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out3}')

print('\nDone! All 3 SVG files saved to figures/')
