import csv, os
base = r'D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\high_freq_analysis'

# Load mapping
mapping_b2a = {}  # bamboo_col → all_col
mapping_a2b = {}  # all_col → bamboo_col
with open(os.path.join(base, 'MSA_Position_Mapping_ALL_to_Bamboo.csv'), 'r', encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        all_pos = int(r['All_Species_MSA_Pos'])
        bm_pos = int(r['Bamboo_MSA_Pos']) if r['Bamboo_MSA_Pos'] != 'None' and r['Bamboo_MSA_Pos'] else None
        if bm_pos:
            mapping_b2a[bm_pos] = all_pos
            mapping_a2b[all_pos] = bm_pos

# L1a sig (ALL_MSA): 25,39,43,44,46,47,69,83,94
# L1b sig (ALL_MSA): 34,43,49,50,83,94
# L2 sig (Bamboo MSA col): 3,26,37,40,41,45,47,48,50,69,70,79,116,118

l1a = {25,39,43,44,46,47,69,83,94}
l1b = {34,43,49,50,83,94}
l2_bamboo = {3,26,37,40,41,45,47,48,50,69,70,79,116,118}

# Map L2 to ALL_MSA
l2_all = set()
for bc in l2_bamboo:
    ac = mapping_b2a.get(bc)
    if ac:
        l2_all.add(ac)
    else:
        print(f'  Bamboo col{bc} → NO ALL_MSA mapping')

print('=== Set sizes ===')
print(f'L1a (All-species): {len(l1a)}')
print(f'L1b (Species-controlled): {len(l1b)}')
print(f'L2 (Bamboo-only mapped): {len(l2_all)} from {len(l2_bamboo)} Bamboo cols')

print(f'\nL1a ∩ L1b: {sorted(l1a & l1b)}')
print(f'L1a ∩ L2: {sorted(l1a & l2_all)}')
print(f'L1b ∩ L2: {sorted(l1b & l2_all)}')
print(f'L1a ∩ L1b ∩ L2: {sorted(l1a & l1b & l2_all)}')
print(f'L1a only: {sorted(l1a - l1b - l2_all)}')
print(f'L1b only: {sorted(l1b - l1a - l2_all)}')
print(f'L2 only: {sorted(l2_all - l1a - l1b)}')

# Full breakdown for 3-set Venn
abc = l1a & l1b & l2_all
ab = (l1a & l1b) - l2_all
ac = (l1a & l2_all) - l1b
bc = (l1b & l2_all) - l1a
a_only = l1a - l1b - l2_all
b_only = l1b - l1a - l2_all
c_only = l2_all - l1a - l1b

print(f'\n=== 3-set Venn values ===')
print(f'A only (L1a): {sorted(a_only)}')
print(f'B only (L1b): {sorted(b_only)}')
print(f'C only (L2): {sorted(c_only)}')
print(f'AB (L1a∩L1b): {sorted(ab)}')
print(f'AC (L1a∩L2): {sorted(ac)}')
print(f'BC (L1b∩L2): {sorted(bc)}')
print(f'ABC (triple): {sorted(abc)}')
