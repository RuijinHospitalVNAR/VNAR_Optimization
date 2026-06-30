"""
Bamboo Shark MSA Fisher检验 — 18条序列 (8 SOLUBLE + 10 INSOLUBLE)
MSA: Bamboo_shark_sequences alignment.1.fa (149列)
排除ENGINEERING，独立分区定义
"""
import pandas as pd, numpy as np
from scipy.stats import fisher_exact
from collections import Counter
import warnings; warnings.filterwarnings('ignore')

msa = ("c:/Users/wcf/.trae-cn/attachments/6a1d3e475a6978b822f2cd3c/"
       "6cbd609c-8718-4aaf-a709-167aa2d1689a_"
       "9009013b-eb26-4b6d-950f-7e9a26198fff_Bamboo_shark_sequences alignment.1.fa")

def load(f): 
    d={}; n=None; c=''
    with open(f) as fh:
        for l in fh:
            l=l.strip()
            if l.startswith('>'): 
                if n: d[n]=c
                n,c=l[1:],''
            elif l: c+=l
        if n: d[n]=c
    return d

all_s = load(msa)
mlen = max(len(s) for s in all_s.values())
print(f"MSA文件: {len(all_s)}条序列, {mlen}列")

# 分类（排除ENGINEERING）
sol_n = sorted([n for n in all_s if 'SOLUBLE' in n.upper() and 'INSOLUBLE' not in n.upper() and 'ENGINEERING' not in n.upper()])
insol_n = sorted([n for n in all_s if 'INSOLUBLE' in n.upper() and 'ENGINEERING' not in n.upper()])
eng_n = sorted([n for n in all_s if 'ENGINEERING' in n.upper()])

sol = [all_s[n] for n in sol_n]
insol = [all_s[n] for n in insol_n]
print(f"\nSOLUBLE: {len(sol)}条")
for n in sol_n: print(f"  {n}")
print(f"INSOLUBLE: {len(insol)}条")
for n in insol_n: print(f"  {n}")
if eng_n: print(f"ENGINEERING(排除): {len(eng_n)}条"); [print(f"  {n}") for n in eng_n]

# ============================================================
# 分区（Bamboo MSA, 149列, 修正版保守残基）
# Cys21/Cys82 = Ig二硫键, Trp35 = CDR1-FR2, Gly112 = CDR3-FR4
# ============================================================
regions = [
    ('FR1',1,21),('CDR1',22,34),('FR2',35,49),('HV2',50,56),
    ('FR3a',57,68),('HV4',69,76),('FR3b',77,82),('CDR3',83,111),('FR4',112,mlen)
]
def region(pos):
    for rn,rs,re in regions:
        if rs<=pos<=re: return rn
    return '?'

# 保守边界验证
print("\n保守边界验证:")
for col,exp,note in [(21,'C','FR1-end'),(35,'W','CDR1-FR2'),(82,'C','FR3b-end'),(112,'G','CDR3-FR4')]:
    aas=[s[col-1] for s in sol+insol if col<=len(s) and s[col-1] not in ('-','.')]
    c=Counter(aas); mc=c.most_common(1)[0]
    ok='✓' if mc[0]==exp and mc[1]/len(aas)>=0.9 else '✗'
    print(f"  {ok} col{col}: {mc[0]} freq={mc[1]/len(aas):.1%} ← {note}")

# ============================================================
# 高频位点（仅Bamboo可溶8条）
# ============================================================
hf = []
for pos in range(mlen):
    aas=Counter(); total=0
    for s in sol:
        if pos<len(s) and s[pos] not in ('-','.'):
            aas[s[pos]]+=1; total+=1
    if not total: continue
    mc,cnt=aas.most_common(1)[0]; freq=cnt/total; cov=total/len(sol)
    if freq>0.70 and cov>0.72:
        hf.append({'position':pos+1,'region':region(pos+1),'aa':mc,'freq':freq,'cov':cov})
print(f"\n高频位点: {len(hf)}个")
for rn,_,_ in regions:
    subset=[x for x in hf if x['region']==rn]
    if subset: print(f"  {rn}: {len(subset)}个 → {[(x['position'],x['aa']) for x in subset[:5]]}{'...' if len(subset)>5 else ''}")

# ============================================================
# Fisher检验（单侧greater, 与原始脚本一致）
# ============================================================
results=[]
for hp in hf:
    pos=hp['position']-1; tgt=hp['aa']; rn=hp['region']
    a=sum(1 for s in sol if len(s)>pos and s[pos]==tgt)
    b=sum(1 for s in sol if len(s)>pos and s[pos] not in ('-','.') and s[pos]!=tgt)
    c=sum(1 for s in insol if len(s)>pos and s[pos]==tgt)
    d=sum(1 for s in insol if len(s)>pos and s[pos] not in ('-','.') and s[pos]!=tgt)
    st,it=a+b,c+d
    if not st or not it: continue
    try:
        _,p=fisher_exact([[a,b],[c,d]],alternative='greater')
        if c==0 and b>0: or_v=float('inf')
        elif a==0 and d>0: or_v=0
        else: or_v=(a*d)/(b*c) if b>0 and c>0 else float('nan')
    except: p,or_v=1.0,float('nan')
    sig=p<0.05
    results.append({'position':pos+1,'region':rn,'target_aa':tgt,
        'sol_match':a,'sol_total':st,'sol_freq':a/st,
        'insol_match':c,'insol_total':it,'insol_freq':c/it,
        'odds_ratio':or_v,'p_value':p,'significant':sig})

df=pd.DataFrame(results)
sig=df[df['significant']]
print(f"\nFisher检验: {len(df)}个位点, 显著: {len(sig)}个")

# ============================================================
# 显著位点展示
# ============================================================
print(f"\n{'='*80}")
print(f"显著位点 (p<0.05): {len(sig)}个")
print(f"{'='*80}")
for _,r in sig.sort_values('p_value').iterrows():
    or_s=f"inf" if r['odds_ratio']==float('inf') else f"{r['odds_ratio']:.1f}"
    print(f"  col{int(r['position']):3d} ({r['target_aa']}) [{r['region']:5s}]: "
          f"sol={r['sol_freq']:.1%}({r['sol_match']}/{r['sol_total']}), "
          f"insol={r['insol_freq']:.1%}({r['insol_match']}/{r['insol_total']}), "
          f"OR={or_s}, p={r['p_value']:.4f}")

# 分区域汇总
print(f"\n分区域:")
for rn,_,_ in regions:
    sub=sig[sig['region']==rn]
    if len(sub):
        print(f"  {rn}: {len(sub)}个 — {list(zip(sub['position'].astype(int),sub['target_aa'].tolist()))}")

# ============================================================
# 仅FR区域的显著位点（与报告对齐）
# ============================================================
fr_sig=sig[sig['region'].str.startswith('FR')]
print(f"\n仅FR区域显著: {len(fr_sig)}个")
for _,r in fr_sig.sort_values('p_value').iterrows():
    or_s=f"inf" if r['odds_ratio']==float('inf') else f"{r['odds_ratio']:.1f}"
    print(f"  Bamboo col{int(r['position']):3d} ({r['target_aa']}) [{r['region']:5s}]: "
          f"sol={r['sol_freq']:.1%}({r['sol_match']}/{r['sol_total']}), "
          f"insol={r['insol_freq']:.1%}({r['insol_match']}/{r['insol_total']}), "
          f"p={r['p_value']:.4f}")

# ============================================================
# 边际显著 (p<0.10)
# ============================================================
mar=df[(df['p_value']<0.10) & (~df['significant'])]
if len(mar):
    print(f"\n边际显著 (p<0.10): {len(mar)}个")
    for _,r in mar.sort_values('p_value').iterrows():
        print(f"  col{int(r['position'])} ({r['target_aa']}) [{r['region']:5s}]: "
              f"sol={r['sol_freq']:.1%}, insol={r['insol_freq']:.1%}, p={r['p_value']:.4f}")

# ============================================================
# 3个robust位点在此MSA中的对应位置（需要映射）
# ============================================================
print(f"\n{'='*80}")
print(f"3个robust位点在全MSA vs Bamboo MSA的对应")
print(f"{'='*80}")
print(f"注意: 全MSA编号≠Bamboo MSA列号, 需通过保守残基映射")
print(f"  全MSA 43(G,FR2) → Bamboo col41(G,FR2)")
print(f"  全MSA 83(I,FR3b) → Bamboo col69(I,HV4区)")
print(f"  全MSA 94(T,FR3b) → Bamboo col79(T,FR3b)")

for label, bcol, exp_aa in [("全MSA-43→Bamboo-41",41,'G'),("全MSA-83→Bamboo-69",69,'I'),("全MSA-94→Bamboo-79",79,'T')]:
    match=df[df['position']==bcol]
    if len(match):
        r=match.iloc[0]
        sig_mark='★' if r['significant'] else ''
        print(f"  {label}({exp_aa}): sol={r['sol_freq']:.1%} insol={r['insol_freq']:.1%} p={r['p_value']:.4f} {sig_mark}")

# 保存
import shutil
out_work = 'c:\\Users\\wcf\\.trae-cn\\work\\6a1d3e475a6978b822f2cd3c\\Bamboo_MSA_18seq_Fisher.csv'
out_dst = 'd:\\文章投递内容\\2026\\VNAR_PRODUCTION_PROTOCOL\\high_freq_analysis\\Bamboo_MSA_18seq_Fisher.csv'
df.to_csv(out_work,index=False)
print(f"\n保存: {out_work}")
try:
    shutil.copy2(out_work, out_dst)
    print(f"已复制到: {out_dst}")
except Exception as e:
    print(f"复制到workspace失败（将使用work目录）: {e}")
print("完成!")
