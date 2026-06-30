"""EC50 FOLD dot plot — Nature palette, SVG output, Illustrator-editable fonts"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ── Nature-inspired palette ──
NATURE_BLUE   = '#3B6FB6'   # round1
NATURE_RED    = '#D62728'   # round2
NATURE_GOLD   = '#EBB434'   # WT diamond
GRID_BLACK    = 'black'
TEXT_DARK     = '#333333'

# ── Fonts: use sans-serif that Illustrator can easily recognize ──
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['mathtext.fontset'] = 'stix'

# ── Data ──
round1 = [
    ('WT',   1.000),
    ('H95A', 0.867),
    ('E96F', 0.534),
    ('S85A', 1.011),
    ('E96Y', 0.591),
    ('P90Y', 0.537),
    ('D88G', 0.619),
]
round2 = [
    ('WT',   1.000),
    ('M01',  2.460),
    ('M02',  1.033),
    ('M07',  0.393),
    ('M10',  0.182),
    ('M05',  0.471),
    ('M08',  0.756),
]

fig, ax = plt.subplots(figsize=(6.5, 6.0))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# ── WT=1 dashed reference line (black, 1pt ≈ 1.3px in matplotlib) ──
ax.axhline(y=1, color='black', linestyle='--', linewidth=1.0, alpha=1.0, zorder=1)

# ── Jittered scatter (bigger dots, no white edge) ──
def plot_group(names_vals, x_pos, color):
    rng = np.random.RandomState(42)
    jitter = rng.uniform(-0.10, 0.10, len(names_vals))
    xs = [x_pos + j for j in jitter]
    ys = [v for _, v in names_vals]

    ax.scatter(xs, ys, c=color, s=130, edgecolors='none',
               zorder=3, alpha=0.90, clip_on=False)

    for i, (name, val) in enumerate(names_vals):
        if val > 1.0:
            ax.annotate(name, (xs[i], val),
                       textcoords="offset points", xytext=(8, 8),
                       fontsize=9, fontweight='bold', color=color,
                       ha='left', va='bottom',
                       arrowprops=dict(arrowstyle='-', color=color, lw=0.6, alpha=0.45))
    return xs, ys

xs1, ys1 = plot_group(round1, 0, NATURE_BLUE)
xs2, ys2 = plot_group(round2, 1, NATURE_RED)

# ── WT diamonds (bigger, no edge) ──
for x, y in [(xs1[0], ys1[0]), (xs2[0], ys2[0])]:
    ax.scatter(x, y, c=NATURE_GOLD, s=160, edgecolors='none',
               zorder=5, marker='D')
ax.text(xs1[0], ys1[0] - 0.055, 'WT', fontsize=8.5, fontweight='bold',
        color=TEXT_DARK, ha='center', va='top', fontfamily='sans-serif')
ax.text(xs2[0], ys2[0] - 0.055, 'WT', fontsize=8.5, fontweight='bold',
        color=TEXT_DARK, ha='center', va='top', fontfamily='sans-serif')

# ── Subtle below-line labels ──
for name, val in round1:
    if name != 'WT' and val < 1.0:
        idx = [n for n, _ in round1].index(name)
        ax.annotate(name, (xs1[idx], val), textcoords="offset points",
                   xytext=(-14, -10), fontsize=6, color=NATURE_BLUE, alpha=0.55,
                   ha='right', va='top')
for name, val in round2:
    if name != 'WT' and val < 1.0:
        idx = [n for n, _ in round2].index(name)
        ax.annotate(name, (xs2[idx], val), textcoords="offset points",
                   xytext=(14, -10), fontsize=6, color=NATURE_RED, alpha=0.55,
                   ha='left', va='top')

# ── Axes ──
ax.set_xlim(-0.30, 1.30)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Round 1\n(Single mutants)', 'Round 2\n(Combinations)'],
                   fontsize=11, fontweight='bold', color=TEXT_DARK, fontfamily='sans-serif')
ax.set_ylabel('EC$_{50}$ (fold over WT)', fontsize=11, color=TEXT_DARK, fontfamily='sans-serif')
ax.set_ylim(0, 2.65)

ax.tick_params(axis='y', labelsize=9, colors=TEXT_DARK)
ax.tick_params(axis='x', pad=10)

# Grid (black, 1pt)
ax.yaxis.grid(True, alpha=1.0, linestyle='-', color='black', linewidth=1.0)
ax.set_axisbelow(True)

# Spine cleanup
for sp in ['top', 'right']:
    ax.spines[sp].set_visible(False)
for sp in ['left', 'bottom']:
    ax.spines[sp].set_color('black')
    ax.spines[sp].set_linewidth(1.0)

plt.tight_layout(pad=1.5)

# ── Save SVG ──
out = r"d:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\high_freq_analysis\EC50_dot_plot.svg"
plt.savefig(out, format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print(f"Saved: {out}")
