import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

fig, ax = plt.subplots(figsize=(7, 6))
v = venn2(
    subsets=(6, 3, 3),
    set_labels=('All-species Fisher\n(9 significant)',
                'Species-controlled\n(6 significant)'),
    set_colors=('#2F5496', '#F4A261'),
    alpha=0.7,
    ax=ax,
)

v.get_label_by_id('10').set_text('6 positions\n25, 39, 44, 46, 47, 69')
v.get_label_by_id('01').set_text('3 positions\n34, 49, 50')
v.get_label_by_id('11').set_text('ROBUST (3)\n43(G)  83(I)  94(T)')

for sid in ['10','01','11']:
    lbl = v.get_label_by_id(sid)
    if lbl:
        lbl.set_fontsize(9)
        lbl.set_fontfamily('Arial')

for sid in ['A','B']:
    lbl = v.get_label_by_id(sid)
    if lbl:
        lbl.set_fontsize(12)
        lbl.set_fontweight('bold')
        lbl.set_fontfamily('Arial')

ax.set_title('Fisher Significance: Cross-Level Overlap\n(ALL_MSA numbering)',
             fontsize=14, fontweight='bold', fontfamily='Arial', pad=20, color='#2F5496')

fig.tight_layout()
out = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\figures\Fig_Chart_Venn_Fisher.svg'
fig.savefig(out, format='svg', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {out}')
