"""Test converter motif extraction fix."""

from app.A_convert_to_hte_format import (
    cached_featurize, 
    format_reaction_key, 
    select_primary_reacted_motifs, 
    select_primary_formed_motif
)
from collections import Counter

# Test Suzuki reaction
smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
reactants, products = smiles.split(">>")
r_list = reactants.split(".")

motif_ids = []
reactant_data = []

print("Testing motif extraction fix...")
print("="*60)

for r in r_list:
    a = cached_featurize(r)
    m = [x.get("id", "") for x in a.get("motifs", []) if x.get("id")]
    reactant_data.append({"motifs": m})
    motif_ids.extend(m)
    print(f"Reactant: {r}")
    print(f"  Motifs: {m}")

p_a = cached_featurize(products)
p_m = [x.get("id", "") for x in p_a.get("motifs", []) if x.get("id")]
print(f"Product: {products}")
print(f"  Motifs: {p_m}")

r_counts = Counter(motif_ids)
p_counts = Counter(p_m)

reacted = {m for m in r_counts if r_counts[m] > p_counts.get(m, 0)}
formed = {m for m in p_counts if p_counts[m] > r_counts.get(m, 0)}
spectators = {m for m in r_counts if r_counts[m] == p_counts.get(m, 0) and r_counts[m] > 0}

print(f"\nReacted: {reacted}")
print(f"Formed: {formed}")
print(f"Spectators: {spectators}")

primary_reacted = select_primary_reacted_motifs(reactant_data, reacted)
primary_formed = select_primary_formed_motif(p_m, formed)

reaction_key = format_reaction_key(
    primary_reacted, 
    [primary_formed] if primary_formed else [], 
    list(spectators)
)

type_a = primary_reacted[0] if len(primary_reacted) > 0 else ""
type_b = primary_reacted[1] if len(primary_reacted) > 1 else ""

print("\n" + "="*60)
print("CSV Column Values:")
print(f"  Reaction_Key: {reaction_key}")
print(f"  reactant_1: {type_a}")
print(f"  reactant_2: {type_b}")

if reaction_key and type_a:
    print("\n✅ FIXED: Columns now have values!")
else:
    print("\n❌ STILL BROKEN: Columns are empty")
