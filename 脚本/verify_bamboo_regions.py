"""
Bamboo MSA 分区边界验证脚本
"""
from collections import Counter

msa_file = (r'c:\Users\wcf\.trae-cn\attachments\6a1d3e475a6978b822f2cd3c'
            r'\6cbd609c-8718-4aaf-a709-167aa2d1689a_'
            r'9009013b-eb26-4b6d-950f-7e9a26198fff_Bamboo_shark_sequences alignment.1.fa')

# Load
d={}; n=None; c=''
with open(msa_file) as f:
    for l in f:
        l=l.strip()
        if l.startswith('>'): 
            if n: d[n]=c
            n,c=l[1:],''
        elif l: c+=l
    if n: d[n]=c

seqs = {n:s for n,s in d.items() if 'ENGINEERING' not in n.upper()}
names = sorted(seqs.keys())
mlen = max(len(s) for s in seqs.values())

regions = [
    ('FR1',1,21),('CDR1',22,34),('FR2',35,49),('HV2',50,56),
    ('FR3a',57,68),('HV4',69,76),('FR3b',77,82),('CDR3',83,111),('FR4',112,mlen)
]
def rgn(p):
    for rn,rs,re in regions:
        if rs<=p<=re: return rn
    return '?'

print("="*90)
print("  Bamboo MSA (18条, 149列) 分区边界验证")
print("="*90)
print(f"  序列数: {len(seqs)} (SOLUBLE=8, INSOLUBLE=10, 排除ENGINEERING=4)")
print()

# ============================================================
# 每列详细分析
# ============================================================
print(f"{'Col':>4s} {'Region':>6s} {'Top AA':>10s} {'Freq':>7s} {'Entropy':>7s}  Notes")
print("-"*70)

boundary_markers = {
    1: '← MSA起始',
    21: '← CONSERVED CYS (Ig disulfide)',
    22: '← CDR1起始',
    35: '← CONSERVED TRP (CDR1→FR2)',
    50: '← HV2起始 (CDR2 equiv)',
    51: '← Gly-Gly motif',
    57: '← FR3a起始',
    69: '← HV4起始 (CDR3-like)',
    77: '← FR3b起始',
    82: '← CONSERVED CYS (Ig disulfide)',
    83: '← CDR3起始',
    84: '← CDR3保守Ala?',
    112: '← CONSERVED GLY (CDR3→FR4)',
    114: '← GGxT motif',
}

for pos in range(mlen):
    col = pos+1
    aas = Counter()
    for s in seqs.values():
        if col <= len(s) and s[col-1] not in ('-','.'):
            aas[s[col-1]] += 1
    total = sum(aas.values())
    if total == 0: continue
    
    mc, cnt = aas.most_common(1)[0]
    freq = cnt/total
    entropy = -sum((c/total)*__import__('math').log2(c/total) for c in aas.values())
    
    note = boundary_markers.get(col, '')
    region_name = rgn(col)
    
    # Show all boundary columns + high-conservation + boundary-adjacent
    if note or freq >= 0.90 or col in [1,21,22,23,34,35,36,49,50,51,56,57,58,67,68,69,70,76,77,78,81,82,83,84,111,112,113,114,115,149]:
        print(f"{col:4d} {region_name:>6s} {mc:>10s} {freq:6.1%} {entropy:7.3f}  {note}")

# ============================================================
# 分区汇总验证
# ============================================================
print()
print("="*90)
print("  分区汇总")
print("="*90)

boundary_checks = [
    (21, 'C', 'FR1→CDR1', 'Ig disulfide Cys'),
    (35, 'W', 'CDR1→FR2', 'Conserved Trp anchor'),
    (82, 'C', 'FR3b→CDR3', 'Ig disulfide Cys (Y-H-C motif)'),
    (112, 'G', 'CDR3→FR4', 'Conserved Gly (GGGT motif start)'),
]

all_ok = True
for col, exp_aa, bound_desc, note in boundary_checks:
    aas = Counter()
    for s in seqs.values():
        if col <= len(s) and s[col-1] not in ('-','.'):
            aas[s[col-1]] += 1
    total = sum(aas.values())
    mc, cnt = aas.most_common(1)[0]
    ok = mc == exp_aa and cnt/total >= 0.90
    
    # Also check the adjacent columns for context
    left_aas = Counter()
    for s in seqs.values():
        if col-1 <= len(s) and s[col-2] not in ('-','.'):
            left_aas[s[col-2]] += 1
    right_aas = Counter()
    for s in seqs.values():
        if col+1 <= len(s) and s[col] not in ('-','.'):
            right_aas[s[col]] += 1
    
    left_top = left_aas.most_common(1)[0] if left_aas else ('?',0)
    right_top = right_aas.most_common(1)[0] if right_aas else ('?',0)
    
    status = '\033[92mOK\033[0m' if ok else '\033[91mFAIL\033[0m'
    print(f"\n  [{status}] col{col} = {mc} ({cnt}/{total} = {cnt/total:.0%})  ← {bound_desc}")
    print(f"       左邻 col{col-1}: {left_top[0]} {left_top[1]}/{sum(left_aas.values())}")
    print(f"       右邻 col{col+1}: {right_top[0]} {right_top[1]}/{sum(right_aas.values())}")
    print(f"       说明: {note}")

# ============================================================
# 区域熵值分析
# ============================================================
print()
print("="*90)
print("  区域保守性统计")
print("="*90)
print(f"{'Region':>6s} {'Cols':>5s} {'Range':>12s} {'Mean Entropy':>13s} {'Top AA cols':>12s}")
print("-"*60)

for rn, rs, re in regions:
    entropies = []
    conserved_count = 0
    for pos in range(rs-1, re):
        aas = Counter()
        for s in seqs.values():
            if pos < len(s) and s[pos] not in ('-','.'):
                aas[s[pos]] += 1
        total = sum(aas.values())
        if total == 0: continue
        freq = max(aas.values())/total
        ent = -sum((c/total)*__import__('math').log2(c/total) for c in aas.values())
        entropies.append(ent)
        if freq >= 0.80:
            conserved_count += 1
    
    mean_ent = sum(entropies)/len(entropies) if entropies else 0
    print(f"{rn:>6s} {re-rs+1:5d} {rs:3d}-{re:<3d}  {mean_ent:13.4f}  {conserved_count:12d}")

# Check Cys at col97
print()
print("="*90)
print("  特别检查: col97 Cys (CDR3内部额外Cys)")
print("="*90)
aas = Counter()
for s in seqs.values():
    if len(s) >= 97 and s[96] not in ('-','.'):
        aas[s[96]] += 1
total = sum(aas.values())
mc, cnt = aas.most_common(1)[0]
print(f"  col97: {mc} {cnt}/{total} ({cnt/total:.0%})")
print(f"  这是CDR3内部的额外Cys，不是Ig二硫键Cys（col82是二硫键Cys）")
print(f"  col97位于CDR3中段（col82+15），保守性高但因为只是Bamboo特征")

print()
print("验证完成!")
