import sys
sys.stdout.reconfigure(encoding='utf-8')

# ============================================================
# L1a positions (ALL_MSA, 9 significant, col34 removed)
# ============================================================
l1a_all = [
    (25, 'FR1', 'R'), (39, 'CDR1', 'R'), (43, 'FR2', 'G'),
    (44, 'FR2', 'S'), (46, 'FR2', 'N'), (47, 'FR2', 'E'),
    (69, 'FR3a', 'V'), (83, 'HV4', 'I'), (94, 'FR3b', 'T'),
]

# L1b positions (ALL_MSA, 6 significant)
l1b_all = [
    (34, 'CDR1', 'S'), (43, 'FR2', 'G'), (49, 'FR2', 'S'),
    (50, 'FR2', 'I'), (83, 'HV4', 'I'), (94, 'FR3b', 'T'),
]

# L2 positions (Bamboo MSA, 14 significant)
l2_bm = [
    (3, 'FR1', 'E'), (26, 'CDR1', 'S'), (37, 'FR2', 'F'),
    (40, 'FR2', 'K'), (41, 'FR2', 'G'), (45, 'FR2', 'K'),
    (47, 'FR2', 'S'), (48, 'FR2', 'L'), (50, 'HV2', 'N'),
    (69, 'HV4', 'I'), (70, 'HV4', 'S'), (79, 'FR3b', 'T'),
    (116, 'FR4', 'I'), (118, 'FR4', 'T'),
]

# ============================================================
# Load ALL→Bamboo position mapping
# ============================================================
import csv, os
mapping_path = r"d:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\high_freq_analysis\MSA_Position_Mapping_ALL_to_Bamboo.csv"
if os.path.exists(mapping_path):
    with open(mapping_path, 'r', encoding='utf-8-sig') as f:
        mapping = list(csv.DictReader(f))
    all_to_bm = {}
    for r in mapping:
        all_to_bm[int(r['All_Species_MSA_Pos'])] = int(r['Bamboo_MSA_Pos'])
    print(f"Loaded {len(all_to_bm)} ALL→Bamboo position mappings")
else:
    all_to_bm = {}
    print("Mapping file not found, using manual mapping")

# ============================================================
# Map L1a and L1b positions to Bamboo MSA
# ============================================================
l1a_bm_set = set()
l1a_bm_detail = []
for pos, region, aa in l1a_all:
    bm_pos = all_to_bm.get(pos)
    if bm_pos:
        l1a_bm_set.add(bm_pos)
        l2_entry = next((r for r in l2_bm if r[0] == bm_pos), None)
        in_l2 = f"YES (L2 col{bm_pos})" if l2_entry else "NO"
        l1a_bm_detail.append((pos, bm_pos, region, aa, in_l2))

l1b_bm_set = set()
l1b_bm_detail = []
for pos, region, aa in l1b_all:
    bm_pos = all_to_bm.get(pos)
    if bm_pos:
        l1b_bm_set.add(bm_pos)
        l2_entry = next((r for r in l2_bm if r[0] == bm_pos), None)
        in_l2 = f"YES (L2 col{bm_pos})" if l2_entry else "NO"
        l1b_bm_detail.append((pos, bm_pos, region, aa, in_l2))

l2_bm_set = set(p for p,_,_ in l2_bm)

# ============================================================
# Compute overlaps
# ============================================================
print(f"\nL1a mapped to Bamboo: {sorted(l1a_bm_set)} ({len(l1a_bm_set)} positions)")
print(f"L1b mapped to Bamboo: {sorted(l1b_bm_set)} ({len(l1b_bm_set)} positions)")
print(f"L2 Bamboo:             {sorted(l2_bm_set)} ({len(l2_bm_set)} positions)")

intersection_all3 = l1a_bm_set & l1b_bm_set & l2_bm_set
l1a_l2_only = (l1a_bm_set & l2_bm_set) - l1b_bm_set
l1b_l2_only = (l1b_bm_set & l2_bm_set) - l1a_bm_set
l1a_l1b_only = (l1a_bm_set & l1b_bm_set) - l2_bm_set
l1a_only = l1a_bm_set - l1b_bm_set - l2_bm_set
l1b_only = l1b_bm_set - l1a_bm_set - l2_bm_set
l2_only = l2_bm_set - l1a_bm_set - l1b_bm_set

print(f"\n--- Overlap Summary (Bamboo MSA numbering) ---")
print(f"L1a ∩ L1b ∩ L2:  {sorted(intersection_all3)} ({len(intersection_all3)} pos)")
print(f"L1a ∩ L2 only:   {sorted(l1a_l2_only)} ({len(l1a_l2_only)} pos)")
print(f"L1b ∩ L2 only:   {sorted(l1b_l2_only)} ({len(l1b_l2_only)} pos)")
print(f"L1a ∩ L1b only:  {sorted(l1a_l1b_only)} ({len(l1a_l1b_only)} pos)")
print(f"L1a only:        {sorted(l1a_only)} ({len(l1a_only)} pos)")
print(f"L1b only:        {sorted(l1b_only)} ({len(l1b_only)} pos)")
print(f"L2 only:         {sorted(l2_only)} ({len(l2_only)} pos)")

# ============================================================
# Detail the triple overlap (robust positions)
# ============================================================
print(f"\n--- Triple overlap (ROBUST) details ---")
for bm_pos in sorted(intersection_all3):
    l1a_entry = next((p for p in l1a_all if all_to_bm.get(p[0]) == bm_pos), None)
    l1b_entry = next((p for p in l1b_all if all_to_bm.get(p[0]) == bm_pos), None)
    l2_entry = next((p for p in l2_bm if p[0] == bm_pos), None)
    l1a_col = l1a_entry[0] if l1a_entry else '?'
    l1b_col = l1b_entry[0] if l1b_entry else '?'
    l2_col = l2_entry[0] if l2_entry else '?'
    l1a_aa = l1a_entry[2] if l1a_entry else '?'
    l2_aa = l2_entry[2] if l2_entry else '?'
    l2_reg = l2_entry[1] if l2_entry else '?'
    print(f"  ALL_MSA col{l1a_col}/{l1b_col} ↔ Bamboo col{bm_pos} ({l2_reg}) AA={l2_aa}")

# ============================================================
# Detail L2-only (bamboo-specific)
# ============================================================
print(f"\n--- L2-only (Bamboo-specific) positions ---")
for bm_pos in sorted(l2_only):
    l2_entry = next((p for p in l2_bm if p[0] == bm_pos), None)
    if l2_entry:
        print(f"  Bamboo col{bm_pos} ({l2_entry[1]}) AA={l2_entry[2]}")

# ============================================================
# Check if col34 is in L1a (it was removed)
# ============================================================
print(f"\n--- Note ---")
print(f"col34(S) in L1b: YES (mapped to Bamboo col{all_to_bm.get(34,'?')})")
print(f"col34(S) in L1a: REMOVED (manually curated out)")
# Check where col34 maps in Bamboo
bm34 = all_to_bm.get(34)
if bm34:
    l2_at_bm34 = next((p for p in l2_bm if p[0] == bm34), None)
    if l2_at_bm34:
        print(f"  col34 ALL → Bamboo col{bm34} → L2 has: col{bm34} ({l2_at_bm34[1]}) AA={l2_at_bm34[2]}")
    else:
        print(f"  col34 ALL → Bamboo col{bm34} → NOT in L2 significant set")
