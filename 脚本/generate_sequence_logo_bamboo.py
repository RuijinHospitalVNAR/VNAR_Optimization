"""
Bamboo MSA Sequence Logo — 基于Bamboo MSA Fisher (14个显著位点, p<0.05)
MSA: Bamboo_shark_sequences alignment.1.fa, 排除ENGINEERING (4条)
分析序列: 8 SOLUBLE + 10 INSOLUBLE = 18条, MSA长度 149列

区域定义 (保守残基分析):
  FR1: 1-21 (Cys21)    CDR1: 22-34 (Trp35)    FR2: 35-49
  HV2: 50-56 (Gly51-52) FR3a: 57-68           HV4: 69-76
  FR3b: 77-82 (Cys82)   CDR3: 83-111          FR4: 112-149 (Gly112)

数据来源: Bamboo_MSA_18seq_Fisher.csv (14个显著位点, 单侧greater, p<0.05)
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt, matplotlib.patches as mpatches
from collections import Counter
import warnings; warnings.filterwarnings('ignore')
plt.rcParams['font.family'] = ['Arial','DejaVu Sans']; plt.rcParams['axes.unicode_minus']=False
plt.rcParams['svg.fonttype']='none'

print("="*70)
print("Bamboo MSA Sequence Logo — 14个Fisher显著位点标记")
print("="*70)

# ============================================================
# MSA文件
# ============================================================
msa_file = ("c:/Users/wcf/.trae-cn/attachments/6a1d3e475a6978b822f2cd3c/"
            "6cbd609c-8718-4aaf-a709-167aa2d1689a_"
            "9009013b-eb26-4b6d-950f-7e9a26198fff_Bamboo_shark_sequences alignment.1.fa")

def load_fasta(fp):
    d={}; n=None; c=''
    with open(fp) as f:
        for l in f:
            l=l.strip()
            if l.startswith('>'): 
                if n: d[n]=c
                n,c=l[1:],''
            elif l: c+=l
        if n: d[n]=c
    return d

all_s=load_fasta(msa_file)
eng=[n for n in all_s if 'ENGINEERING' in n.upper()]
filt={n:s for n,s in all_s.items() if 'ENGINEERING' not in n.upper()}
seqs=list(filt.values()); n_seq=len(seqs); msa_len=max(len(s) for s in seqs)
sol_n=sorted([n for n in filt if 'SOLUBLE' in n.upper() and 'INSOLUBLE' not in n.upper()])
insol_n=sorted([n for n in filt if 'INSOLUBLE' in n.upper()])
print(f"MSA: {len(all_s)}条 → 排除ENGINEERING {len(eng)}条 → {n_seq}条 (SOL={len(sol_n)}, INSOL={len(insol_n)}), {msa_len}列")

# ============================================================
# 读取Bamboo MSA Fisher结果（14个显著位点, 1-sided greater）
# ============================================================
fisher_csv_work = ('c:\\Users\\wcf\\.trae-cn\\work\\6a1d3e475a6978b822f2cd3c\\'
                   'Bamboo_MSA_18seq_Fisher.csv')
fisher_csv_ws = ('d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\'
                 'high_freq_analysis\\Bamboo_MSA_18seq_Fisher.csv')
import os as _os
fisher_csv = fisher_csv_work if _os.path.exists(fisher_csv_work) else fisher_csv_ws
fisher_df = pd.read_csv(fisher_csv)
sig = fisher_df[fisher_df['significant']]
sig_positions = set(zip(sig['position'].astype(int), sig['target_aa']))
sig_cols = set(sig['position'].astype(int))
print(f"Bamboo Fisher显著位点: {len(sig)}个")
for _,r in sig.sort_values('p_value').iterrows():
    print(f"  col{int(r['position']):3d} ({r['target_aa']}) [{r['region']}]: "
          f"sol={r['sol_freq']:.1%}, insol={r['insol_freq']:.1%}, p={r['p_value']:.4f}")

# ============================================================
# 区域定义
# ============================================================
regions=[
    ('FR1',1,21,'#E8E8E8','#333'),('CDR1',22,34,'#FFE4E4','#C00'),
    ('FR2',35,49,'#E8E8E8','#333'),('HV2',50,56,'#E8D4F0','#93C'),
    ('FR3a',57,68,'#E8E8E8','#333'),('HV4',69,76,'#E8D4F0','#93C'),
    ('FR3b',77,82,'#E8E8E8','#333'),('CDR3',83,111,'#FFE4E4','#C00'),
    ('FR4',112,msa_len,'#E8E8E8','#333')
]
def rgn(pos):
    for rn,rs,re,*_ in regions:
        if rs<=pos<=re: return rn
    return '?'

print("\n区域定义:")
for rn,rs,re,*_ in regions:
    tag=""
    if rn=='FR1': tag=" ← Cys21"
    elif rn=='CDR1': tag=" ← Trp35"
    elif rn=='FR2': tag=" ← Trp35+"
    elif rn=='FR3b': tag=" ← Cys82"
    elif rn=='CDR3': tag=" ← C82+1→G111"
    elif rn=='FR4': tag=" ← Gly112"
    print(f"  {rn:6s}: {rs:3d}-{re:3d} ({re-rs+1:3d}){tag}")

# ============================================================
# 氨基酸颜色
# ============================================================
aa_colors = {
    'G':'#CCCCCC','P':'#FF69B4',
    'A':'#2D2D2D','V':'#2D2D2D','L':'#2D2D2D','I':'#2D2D2D','M':'#2D2D2D',
    'F':'#2D2D2D','W':'#2D2D2D','Y':'#2D2D2D','C':'#FFD700',
    'S':'#00A86B','T':'#00A86B','N':'#00A86B','Q':'#00A86B',
    'K':'#1E90FF','R':'#1E90FF','H':'#1E90FF',
    'D':'#DC143C','E':'#DC143C',
}

def calc_info(counts, n_total):
    tot=sum(counts.values())
    if not tot: return 0
    probs=[c/tot for c in counts.values() if c>0]
    ent=-sum(p*np.log2(p) for p in probs)
    return (np.log2(20)-ent)*(tot/n_total)

# ============================================================
# 计算氨基酸分布 + 保守边界验证
# ============================================================
print("\n保守边界验证:")
for col,exp,note in [(21,'C','FR1-end'),(35,'W','CDR1-FR2'),(82,'C','FR3b-end'),(112,'G','CDR3-FR4')]:
    aas=[s[col-1] for s in seqs if col<=len(s) and s[col-1] not in ('-','.')]
    c=Counter(aas); mc=c.most_common(1)[0]
    ok='✓' if mc[0]==exp and mc[1]/len(aas)>=.9 else '✗'
    print(f"  {ok} col{col}: {mc[0]} freq={mc[1]/len(aas):.1%} ← {note}")

pos_data=[]
for pos in range(msa_len):
    aa_counts=Counter()
    for s in seqs:
        if pos<len(s) and s[pos]!='-': aa_counts[s[pos]]+=1
    tot=sum(aa_counts.values())
    info=calc_info(aa_counts, n_seq)
    pos_data.append({'position':pos+1,'counts':aa_counts,'total':tot,'info':info,
                     'is_sig':(pos+1) in sig_cols})

# ============================================================
# Sequence Logo 绘图
# ============================================================
print("生成 Sequence Logo...")
# A4竖版宽度: 210mm ≈ 8.27 inches, 高度减半
fig,ax=plt.subplots(figsize=(8.27,2.25))
bar_width=.85

for i,d in enumerate(pos_data[:msa_len]):
    if not d['total']: continue
    sorted_aa=sorted(d['counts'].items(),key=lambda x:x[1])
    y_bottom=0
    for aa,cnt in sorted_aa:
        freq=cnt/d['total']; height=freq*d['info']
        color=aa_colors.get(aa,'#808080')
        alpha=1.0 if d['is_sig'] else 0.9
        ax.bar(i+.5,height,bottom=y_bottom,width=bar_width,color=color,edgecolor='none',alpha=alpha)
        if freq>0.5:
            tc='white' if color in ['#2D2D2D','#1E90FF','#DC143C','#FFD700'] else 'black'
            ax.text(i+.5,y_bottom+height/2,aa,ha='center',va='center',fontsize=6,
                   fontweight='bold',color=tc)
        y_bottom+=height

ax.set_xlim(0,msa_len); ax.set_ylim(0,4.5)
ax.set_ylabel('Information (bits)',fontsize=9,fontweight='bold')
ax.set_title('VNAR Sequence Logo — Bamboo Shark MSA (SOLUBLE+INSOLUBLE, n=18, 149 columns)\nFR1-CDR1-FR2-HV2-FR3a-HV4-FR3b-CDR3-FR4',fontsize=11,fontweight='bold',pad=10)

# 位置编号（每5列标注）
for p in range(1,msa_len+1,5):
    ax.plot([p-.5,p-.5],[0,-.06],color='#999',lw=.3)
    ax.text(p-.5,-.20,str(p),rotation=90,ha='center',va='top',fontsize=8,color='#333')
ax.plot([1-.5,1-.5],[0,-.06],color='#999',lw=.3)
ax.text(1-.5,-.20,'1',rotation=90,ha='center',va='top',fontsize=8,color='#333')

# 区域标注
y_line=-.70; lh=.15
for i,(rn,rs,re,bg,lc) in enumerate(regions):
    if rs>msa_len: break
    re2=min(re,msa_len)
    ax.axvspan(rs-1.5,re2-.5,alpha=.25,color=bg)
    ax.plot([rs-1.5,re2-.5],[y_line,y_line],color=lc,lw=1.5)
    if i>0:
        p2=min(regions[i-1][2],msa_len)
        ax.plot([p2-.5,p2-.5],[y_line-lh,y_line+lh],color='#666',lw=.6)
    xc=(rs-1.5+re2-.5)/2; w=re2-rs+1
    fs=7 if w>20 else 6 if w>10 else 5 if w>5 else 4
    ax.text(xc,y_line-.28,rn,ha='center',va='top',fontsize=fs,fontweight='bold',color=lc)

# Fisher显著位点标记（▼+AA标注）
y_sig=4.35
for _,r in sig.iterrows():
    col=int(r['position']); tgt_aa=r['target_aa']
    if col<=msa_len:
        ax.plot(col-.5,y_sig,'v',markersize=4,color='#FF6600',alpha=.95,mew=.6,mec='#CC3300')
        ax.text(col-.5,y_sig+.10,tgt_aa,ha='center',va='bottom',fontsize=4,
               color='#CC3300',fontweight='bold')

# 保守Cys标记（★）
for p in [21,82]:
    if p<=msa_len:
        ax.plot(p-.5,4.05,'*',markersize=7,color='#FFD700',mew=.6,mec='#B8860B')

ax.set_ylim(-1.0,4.6); ax.set_xticks([])

# 图例
leg1=[mpatches.Patch(facecolor='#E8E8E8',edgecolor='#999',label='FR (Framework)',alpha=.5),
      mpatches.Patch(facecolor='#FFE4E4',edgecolor='#C00',label='CDR',alpha=.5),
      mpatches.Patch(facecolor='#E8D4F0',edgecolor='#93C',label='HV (Hypervariable)',alpha=.5)]
leg2=[mpatches.Patch(facecolor='#2D2D2D',label='Hydrophobic'),
      mpatches.Patch(facecolor='#00A86B',label='Polar'),
      mpatches.Patch(facecolor='#1E90FF',label='Basic'),
      mpatches.Patch(facecolor='#DC143C',label='Acidic')]
leg3=[plt.Line2D([0],[0],marker='v',color='w',markerfacecolor='#FF6600',markersize=8,
        label='Fisher signif. (p<0.05)'),
      plt.Line2D([0],[0],marker='*',color='w',markerfacecolor='#FFD700',markersize=10,
        label='Conserved Cys (Ig disulfide)')]
fig.legend(handles=leg1+leg2+leg3,loc='upper right',fontsize=6,title='Legend',title_fontsize=7,
          frameon=True,bbox_to_anchor=(0.98,0.98))
plt.tight_layout()

# 保存
out_svg='d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\figures\\VNAR_SequenceLogo_Bamboo_149col.svg'
out_pdf='d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\figures\\VNAR_SequenceLogo_Bamboo_149col.pdf'
plt.savefig(out_svg,format='svg',dpi=300); plt.savefig(out_pdf,format='pdf',dpi=300)
plt.close()
print(f"\nSVG: {out_svg}")
print(f"PDF: {out_pdf}")
print("完成!")

# 将脚本复制到high_freq_analysis
import shutil as _sh
dst='c:\\Users\\wcf\\.trae-cn\\work\\6a1d3e475a6978b822f2cd3c\\script_sequence_logo_bamboo.py'
_sh.copy2(__file__,dst)
print(f"脚本副本: {dst}")
try:
    dst_ws='d:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\high_freq_analysis\\script_sequence_logo_bamboo.py'
    _sh.copy2(__file__,dst_ws)
    print(f"脚本已复制到workspace: {dst_ws}")
except Exception as e:
    print(f"复制到workspace失败: {e}")
