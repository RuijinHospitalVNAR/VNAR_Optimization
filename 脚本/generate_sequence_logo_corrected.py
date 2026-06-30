"""
生成VNAR Sequence Logo图 — 基于MSA保守残基分析修正区域边界（v2）
原Feng et al. 2019的130位编号系统问题已修正。
修正后区域定义（基于260-sequence MSA保守残基分析，205列）。
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.family'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'

print("=" * 70)
print("VNAR Sequence Logo — MSA保守残基修正版 (全205列)")
print("=" * 70)

# ============================================================
# MSA文件
# ============================================================
msa_file = ('c:\\Users\\wcf\\.trae-cn\\attachments\\6a1d3e475a6978b822f2cd3c\\'
            '75f0990d-adc6-412b-8348-059a060267fd_'
            '68d83b87-7b75-4b53-b501-71d24f30e076_VNAR_engneer_for_MSA alignment.fa')

def read_multiline_fasta(fp):
    seqs, name, seq = {}, None, ''
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name: seqs[name] = seq
                name, seq = line[1:], ''
            elif line: seq += line
        if name: seqs[name] = seq
    return seqs

all_seqs = read_multiline_fasta(msa_file)
sequences = list(all_seqs.values())
alignment_length = len(sequences[0])
print(f"序列数: {len(sequences)}, MSA长度: {alignment_length}")

# ============================================================
# 辅助数据
# ============================================================
hf_file = ('d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\'
           'high_freq_analysis\\MSA_high_frequency_positions_final_corrected.csv')
hf_df = pd.read_csv(hf_file)
fisher_file = ('d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\'
               'high_freq_analysis\\FR_fisher_test_results.csv')
fisher_df = pd.read_csv(fisher_file)
sig_positions = fisher_df[fisher_df['significant']]['position'].tolist()
print(f"高频位点: {len(hf_df)}, Fisher显著: {len(sig_positions)}")

# ============================================================
# 修正后的VNAR区域定义（MSA列号，基于保守残基分析）
#   Cys22 (99.2%) → FR1末尾
#   Trp37 (99.2%) → CDR1-FR2边界
#   GG 59-60     → HV2-FR3a边界
#   Cys98 (98.5%) → FR3b-CDR3边界
#   YGGGT 160-164 → CDR3-FR4边界 (WGxGT-like)
# ============================================================
regions = [
    ('FR1',  1,  22, '#E8E8E8', '#333333'),
    ('CDR1', 23, 36, '#FFE4E4', '#CC0000'),
    ('FR2',  37, 50, '#E8E8E8', '#333333'),
    ('HV2',  51, 58, '#E8D4F0', '#9932CC'),
    ('FR3a', 59, 70, '#E8E8E8', '#333333'),
    ('HV4',  71, 78, '#E8D4F0', '#9932CC'),
    ('FR3b', 79, 98, '#E8E8E8', '#333333'),
    ('CDR3', 99, 159, '#FFE4E4', '#CC0000'),
    ('FR4',  160, 205, '#E8E8E8', '#333333'),
]

print("\n修正后区域定义（MSA列号）:")
for name, start, end, _, _ in regions:
    tag = ""
    if name == 'FR1': tag = " ← Cys22 (99.2%)"
    elif name == 'CDR1': tag = " ← Trp37边界"
    elif name == 'FR2': tag = " ← Trp37开始"
    elif name == 'FR3b': tag = " ← Cys98 (98.5%)"
    elif name == 'FR4': tag = " ← YGGGT 160-164"
    print(f"  {name:6s}: {start:3d}-{end:3d} ({end-start+1:3d}){tag}")

# ============================================================
# 氨基酸颜色
# ============================================================
aa_colors = {
    'G': '#CCCCCC', 'P': '#FF69B4',
    'A': '#2D2D2D', 'V': '#2D2D2D', 'L': '#2D2D2D', 'I': '#2D2D2D', 'M': '#2D2D2D',
    'F': '#2D2D2D', 'W': '#2D2D2D', 'Y': '#2D2D2D', 'C': '#FFD700',
    'S': '#00A86B', 'T': '#00A86B', 'N': '#00A86B', 'Q': '#00A86B',
    'K': '#1E90FF', 'R': '#1E90FF', 'H': '#1E90FF',
    'D': '#DC143C', 'E': '#DC143C',
}

def calc_info_content(counts, total_seqs):
    """计算信息含量：log2(20) - entropy，再乘以覆盖率（避免gap污染）"""
    total = sum(counts.values())
    if total == 0: return 0
    probs = [c / total for c in counts.values() if c > 0]
    entropy = -sum(p * np.log2(p) for p in probs)
    coverage = total / total_seqs
    return (np.log2(20) - entropy) * coverage

# ============================================================
# 计算每个位置的氨基酸分布
# ============================================================
print("\n计算氨基酸分布...")
position_data = []
for pos in range(alignment_length):
    aa_counts = Counter()
    for seq in sequences:
        if pos < len(seq) and seq[pos] != '-':
            aa_counts[seq[pos]] += 1
    total = sum(aa_counts.values())
    info = calc_info_content(aa_counts, len(sequences))
    hf_match = hf_df[hf_df['position'] == pos + 1]
    position_data.append({
        'position': pos + 1,
        'counts': aa_counts, 'total': total,
        'info_content': info,
        'is_hf': len(hf_match) > 0,
        'hf_aa': hf_match.iloc[0]['most_common_aa'] if len(hf_match) > 0 else None,
        'is_sig': pos + 1 in sig_positions,
    })

# ============================================================
# Sequence Logo 绘图
# ============================================================
print("生成 Sequence Logo...")
n_positions = alignment_length  # 205列，全图

# A4竖版宽度: 210mm ≈ 8.27 inches, 高度减半
fig_width = 8.27
fig_height = 2.25
fig = plt.figure(figsize=(fig_width, fig_height))
ax = fig.add_subplot(111)

bar_width = 0.85

for i, data in enumerate(position_data[:n_positions]):
    counts = data['counts']
    total = data['total']
    info = data['info_content']
    if total == 0: continue
    sorted_aa = sorted(counts.items(), key=lambda x: x[1])
    y_bottom = 0
    for aa, count in sorted_aa:
        freq = count / total
        height = freq * info
        color = aa_colors.get(aa, '#808080')
        alpha = 1.0 if data['is_sig'] else 0.9
        ax.bar(i + 0.5, height, bottom=y_bottom, width=bar_width,
               color=color, edgecolor='none', alpha=alpha)
        if freq > 0.5:
            txt_color = 'white' if color in ['#2D2D2D', '#1E90FF', '#DC143C', '#FFD700'] else 'black'
            ax.text(i + 0.5, y_bottom + height / 2, aa,
                   ha='center', va='center', fontsize=6,
                   fontweight='bold', color=txt_color)
        y_bottom += height

ax.set_xlim(0, n_positions)
ax.set_ylim(0, 4.5)
ax.set_ylabel('Information (bits)', fontsize=9, fontweight='bold')
ax.set_title('VNAR Sequence Logo — 260-Sequence MSA (205 columns)\nFR1-CDR1-FR2-HV2-FR3a-HV4-FR3b-CDR3-FR4', 
            fontsize=11, fontweight='bold', pad=10)

# ============================================================
# 位置编号（每5列标注）
# ============================================================
def add_position_numbers(ax, n, interval=5):
    y_num = -0.20
    for p in range(1, n + 1):
        if (p - 1) % interval == 0 or p == 1:
            ax.plot([p - 0.5, p - 0.5], [0, -0.06], color='#999', linewidth=0.3)
            ax.text(p - 0.5, y_num, str(p), rotation=90, ha='center', va='top',
                   fontsize=8, color='#333')

# ============================================================
# 区域标注
# ============================================================
def add_region_labels(ax, n, regions):
    y_line = -0.70
    line_h = 0.15
    for i, (name, rs, re, bg, lc) in enumerate(regions):
        if rs > n: break
        re2 = min(re, n)
        ax.axvspan(rs - 1.5, re2 - 0.5, alpha=0.25, color=bg)
        ax.plot([rs - 1.5, re2 - 0.5], [y_line, y_line], color=lc, linewidth=1.5)
        if i > 0:
            prev = regions[i-1]
            ax.plot([prev[1] - 1.5, prev[1] - 1.5],
                   [y_line - line_h, y_line + line_h], color='#666', linewidth=0.6)
        label_x = (rs - 1.5 + re2 - 0.5) / 2
        w = re2 - rs + 1
        fs = 7 if w > 20 else 6 if w > 10 else 5 if w > 5 else 4
        ax.text(label_x, y_line - 0.28, name, ha='center', va='top',
               fontsize=fs, fontweight='bold', color=lc)

# ============================================================
# 标记
# ============================================================
hf_positions = [p for p in hf_df['position'].tolist() if 1 <= p <= n_positions]
y_mark_top = 4.35
y_cys_top = 4.05

for p in hf_positions:
    ax.plot(p - 0.5, y_mark_top, 'v', markersize=2.5, color='#FF6600', alpha=0.9)

for p in [22, 98]:
    if p <= n_positions:
        ax.plot(p - 0.5, y_cys_top, '*', markersize=5, color='#FFD700', alpha=1.0)

# ============================================================
# 应用标注
# ============================================================
add_position_numbers(ax, n_positions, interval=5)
add_region_labels(ax, n_positions, regions)

ax.set_ylim(-1.0, 4.5)
ax.set_xticks([])

# ============================================================
# 图例
# ============================================================
legend_elements = [
    mpatches.Patch(facecolor='#E8E8E8', edgecolor='#333', label='FR', alpha=0.5),
    mpatches.Patch(facecolor='#FFE4E4', edgecolor='#C00', label='CDR', alpha=0.5),
    mpatches.Patch(facecolor='#E8D4F0', edgecolor='#93C', label='HV', alpha=0.5),
    mpatches.Patch(facecolor='#FFD700', edgecolor='#B8860B', label='Cys', alpha=0.8),
]
aa_legend = [
    mpatches.Patch(facecolor='#2D2D2D', label='Hydrophobic'),
    mpatches.Patch(facecolor='#00A86B', label='Polar'),
    mpatches.Patch(facecolor='#1E90FF', label='Basic'),
    mpatches.Patch(facecolor='#DC143C', label='Acidic'),
]
marker_legend = [
    plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='#FF6600', markersize=8, label='High freq.'),
    plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='#FFD700', markersize=10, label='Conserved Cys'),
]

fig.legend(handles=legend_elements + aa_legend + marker_legend, loc='upper right',
          fontsize=6, title='Legend', title_fontsize=7, frameon=True,
          bbox_to_anchor=(0.98, 0.98))

plt.tight_layout()

# ============================================================
# 保存
# ============================================================
output_svg = ('d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\'
              'figures\\VNAR_SequenceLogo_205col_Corrected.svg')
output_pdf = ('d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\'
              'figures\\VNAR_SequenceLogo_205col_Corrected.pdf')

plt.savefig(output_svg, format='svg', dpi=300)
plt.savefig(output_pdf, format='pdf', dpi=300)
plt.close()

print(f"SVG: {output_svg}")
print(f"PDF: {output_pdf}")
print("\n" + "=" * 70)
print("完成！区域边界已基于MSA保守残基修正")
print("=" * 70)
