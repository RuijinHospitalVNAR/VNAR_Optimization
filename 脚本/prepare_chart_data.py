import csv, os, json
from collections import Counter

base = r"D:\文章投递内容\2026\VNAR_PRODUCTION_PROTOCOL\high_freq_analysis"

# ============================================================
# Chart 1: Species + Label distribution
# ============================================================

# Read species data from JSON
with open(os.path.join(base, "species_analysis_complete.json"), "r", encoding="utf-8") as f:
    jdata = json.load(f)

species_counts = {}
for sp in jdata.get("species_order", []):
    spd = jdata.get("species_data", {}).get(sp, {})
    species_counts[sp] = len(spd.get("sequences", []))

unannotated = 256 - sum(species_counts.values())

# Read MSA to count SOLUBLE/INSOLUBLE
fn_msa = os.path.join(base, "VNAR_engneer_for_MSA alignment.fa")
d = {}
n = None
with open(fn_msa) as f:
    for l in f:
        l = l.strip()
        if l.startswith('>'):
            if n: d[n] = c
            n, c = l[1:], ''
        elif l: c += l
    if n: d[n] = c

sol_count = 0; ins_count = 0; unl_count = 0
for name in d:
    if 'ENGINEERING' in name.upper():
        continue
    full = name.upper()
    if 'SOLUBLE' in full:
        sol_count += 1
    elif 'INSOLUBLE' in full:
        ins_count += 1
    else:
        unl_count += 1

# ============================================================
# Chart 1 Data: Multi-level donut
# ============================================================
print("=== CHART 1: Species & Label Distribution ===")
print("Species (outer ring):")
for sp, cnt in sorted(species_counts.items(), key=lambda x: -x[1]):
    print(f"  {sp}: {cnt} ({cnt/256*100:.1f}%)")
print(f"  Unannotated: {unannotated} ({unannotated/256*100:.1f}%)")

print(f"\nLabel status:")
print(f"  SOLUBLE: {sol_count}")
print(f"  INSOLUBLE: {ins_count}")
print(f"  UNLABELED: {unl_count}")

# ============================================================
# Chart 2: Venn - 3 analysis levels overlap
# ============================================================
# Level 1a sig (ALL_MSA): 25,39,43,44,46,47,69,83,94
# Level 1b sig (ALL_MSA): 34,43,49,50,83,94
# Level 2 sig (Bamboo MSA): different numbering, need mapping

# Load mapping
mapping = {}
with open(os.path.join(base, "MSA_Position_Mapping_ALL_to_Bamboo.csv"), "r", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        all_pos = int(r["All_Species_MSA_Pos"])
        bm_pos = int(r["Bamboo_MSA_Pos"]) if r["Bamboo_MSA_Pos"] != "None" else None
        mapping[all_pos] = bm_pos

l1a = {25,39,43,44,46,47,69,83,94}
l1b = {34,43,49,50,83,94}
l2_bamboo_cols = {3,26,37,40,41,45,47,48,50,69,70,79,116,118}

# Map Level 2 to ALL_MSA
l2_mapped = set()
for bm_col in l2_bamboo_cols:
    for all_pos, bm_pos in mapping.items():
        if bm_pos == bm_col:
            l2_mapped.add(all_pos)

print(f"\n=== CHART 2: Venn Diagram Data ===")
print(f"Level 1a (All-species): {len(l1a)} sig positions")
print(f"Level 1b (Species-controlled): {len(l1b)} sig positions")
print(f"Level 2 (Bamboo-only, mapped to ALL_MSA): {len(l2_mapped)} sig positions")
print(f"\nOverlaps:")
print(f"  L1a ∩ L1b: {l1a & l1b}")
print(f"  L1a ∩ L2: {l1a & l2_mapped}")
print(f"  L1b ∩ L2: {l1b & l2_mapped}")
print(f"  L1a ∩ L1b ∩ L2: {l1a & l1b & l2_mapped}")

# Save as JSON for chart generation
chart1_data = {
    "species": [{"category": k, "value": v} for k, v in sorted(species_counts.items(), key=lambda x: -x[1])] + 
               [{"category": "Unannotated", "value": unannotated}],
    "labels": [{"category": "SOLUBLE (" + str(sol_count) + ")", "value": sol_count},
               {"category": "INSOLUBLE (" + str(ins_count) + ")", "value": ins_count},
               {"category": "UNLABELED (" + str(unl_count) + ")", "value": unl_count}]
}

# For Venn: use Level 1a and Level 1b overlap (same numbering)
overlap = l1a & l1b
l1a_only = l1a - l1b
l1b_only = l1b - l1a

chart2_data = {
    "l1a": {
        "set_name": "All-species Fisher",
        "positions": sorted(list(l1a)),
        "count": len(l1a)
    },
    "l1b": {
        "set_name": "Species-controlled",
        "positions": sorted(list(l1b)),
        "count": len(l1b)
    },
    "overlap": {
        "positions": sorted(list(overlap)),
        "count": len(overlap)
    }
}

with open(os.path.join(base, "chart_data.json"), "w") as f:
    json.dump({"chart1": chart1_data, "chart2": chart2_data}, f, indent=2)

print(f"\nChart data saved to chart_data.json")
print(f"\n=== For Chart 2 Venn (2-set, ALL_MSA numbering) ===")
print(f"  L1a only: {sorted(l1a_only)}")
print(f"  L1b only: {sorted(l1b_only)}")
print(f"  Overlap: {sorted(overlap)}")
