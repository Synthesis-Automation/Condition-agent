"""Debug reaction key extraction."""

from chemtools.featurizers.formatters.reaction import featurize_reaction

# Test the dithiocarbamate reaction
rxn = "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1"

result = featurize_reaction(rxn)
aggregates = result.get("aggregates", {})

print("=== Reaction Key Debug ===\n")
print(f"Reaction: {rxn[:50]}...\n")

print("Aggregates:")
print(f"  Reacted: {aggregates.get('reacted_motifs', [])}")
print(f"  Formed: {aggregates.get('formed_motifs', [])}")
print(f"  Spectator: {aggregates.get('spectator_motifs', [])}")

print(f"\nReaction Key: {result.get('reaction_key')}")

# Now test the matcher directly
from chemtools.reaction_key_matcher import detect_from_reaction_key, _load_compound_logic

reacted = aggregates.get('reacted_motifs', [])
formed = aggregates.get('formed_motifs', [])

print(f"\nCalling detect_from_reaction_key with:")
print(f"  reacted: {reacted}")
print(f"  formed: {formed}")

top, matches = detect_from_reaction_key(reacted, formed)

if top:
    print(f"\nTop match: {top.reaction_type} (conf: {top.confidence})")
else:
    print("\nNo matches found!")
    
# Debug: Check what sets are loaded
logic = _load_compound_logic()
motif_sets = logic.get("motif_sets", {})
print(f"\nLoaded motif sets: {list(motif_sets.keys())}")

# Check sp2_electrophiles
sp2 = motif_sets.get("sp2_electrophiles", {}).get("members", [])
print(f"\nsp2_electrophiles sample: {sp2[:5]}...")

# Check if Ar-I is in sp2_electrophiles
if "Ar-I" in sp2:
    print("✓ Ar-I is in sp2_electrophiles")
else:
    print("✗ Ar-I NOT in sp2_electrophiles")

# Check thiols_sh
thiols = motif_sets.get("thiols_sh", {}).get("members", [])
print(f"\nthiols_sh: {thiols}")

# Check if Thiocarbonyl-SH is in thiols_sh
if "Thiocarbonyl-SH" in thiols:
    print("✓ Thiocarbonyl-SH is in thiols_sh")
else:
    print("✗ Thiocarbonyl-SH NOT in thiols_sh")
