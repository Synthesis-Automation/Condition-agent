"""Debug why Suzuki detection is failing for the given reaction."""

from chemtools.featurizers.detection import detect_reaction_types
from chemtools.featurizers.detection.matchers import _detect_motif_profile

# Your reaction
reaction_smiles = "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1"

print("=" * 80)
print("REACTION SMILES:", reaction_smiles)
print("=" * 80)

# Parse reactants and products
parts = reaction_smiles.split(">>")
reactants_str = parts[0]
products_str = parts[1] if len(parts) > 1 else ""

reactants = reactants_str.split(".")
products = products_str.split(".") if products_str else []

print("\nREACTANTS:")
for i, r in enumerate(reactants):
    print(f"  {i}: {r}")

print("\nPRODUCTS:")
for i, p in enumerate(products):
    print(f"  {i}: {p}")

print("\n" + "=" * 80)
print("MOTIF DETECTION IN REACTANTS")
print("=" * 80)

# Detect motifs in reactants
reactant_profile = _detect_motif_profile(reactants)
print(f"\nDetected {len(reactant_profile)} unique motifs in reactants:")
for motif_id, data in sorted(reactant_profile.items()):
    print(f"  - {motif_id}: count={data['count']}, molecules={sorted(data['molecules'])}")

print("\n" + "=" * 80)
print("MOTIF DETECTION IN PRODUCTS")
print("=" * 80)

# Detect motifs in products
product_profile = _detect_motif_profile(products)
print(f"\nDetected {len(product_profile)} unique motifs in products:")
for motif_id, data in sorted(product_profile.items()):
    print(f"  - {motif_id}: count={data['count']}, molecules={sorted(data['molecules'])}")

print("\n" + "=" * 80)
print("REACTION TYPE DETECTION")
print("=" * 80)

# Detect reaction types (without product confirmation)
result = detect_reaction_types(reaction_smiles, confirm_coupling_products=False)

print(f"\nFound {len(result.matches)} match(es):")
for i, match in enumerate(result.matches):
    print(f"\n  Match {i + 1}:")
    print(f"    Type: {match.reaction_type}")
    print(f"    Name: {match.name}")
    print(f"    Confidence: {match.matched_slots}/{match.required_slots}")
    print(f"    Slot Evidence:")
    for slot, motifs in match.slot_evidence.items():
        print(f"      {slot}: {motifs}")

print("\n" + "=" * 80)
print("SUZUKI REQUIREMENTS (from reaction_types.v4.0.json)")
print("=" * 80)
print("""
Suzuki_miyaura:
  reactants:
    electrophile: @sp2_electrophiles
      - Ar-I, Ar-Br, Ar-Cl, Ar-F
      - Ar-OTf, Ar-OMs, Ar-OTs
      - Alkenyl-I, Alkenyl-Br, etc.
    nucleophile: @organoboron
      - Ar-B(OH)2, Ar-Bpin, Ar-BF3K, Ar-B(OR)2
      - Alkenyl-B(OH)2, Alkenyl-Bpin, etc.
  products:
    product: [Ar-Ar, Ar-Alkenyl, Alkenyl-Alkenyl]
""")

print("\n" + "=" * 80)
print("ANALYSIS")
print("=" * 80)

# Check which Suzuki motifs are present
sp2_electrophiles = [
    "Ar-I", "Ar-Br", "Ar-Cl", "Ar-F",
    "Ar-OTf", "Ar-OMs", "Ar-OTs",
    "Alkenyl-I", "Alkenyl-Br", "Alkenyl-Cl", "Alkenyl-F",
    "Alkenyl-OTf", "Alkenyl-OMs", "Alkenyl-OTs"
]

organoboron = [
    "Ar-B(OH)2", "Ar-Bpin", "Ar-BF3K", "Ar-B(OR)2",
    "Alkenyl-B(OH)2", "Alkenyl-Bpin", "Alkenyl-BF3K", "Alkenyl-B(OR)2"
]

coupling_products = ["Ar-Ar", "Ar-Alkenyl", "Alkenyl-Alkenyl"]

print("\nChecking for sp2_electrophiles:")
found_electrophiles = [m for m in sp2_electrophiles if m in reactant_profile]
if found_electrophiles:
    print(f"  ✓ FOUND: {found_electrophiles}")
else:
    print("  ✗ NOT FOUND")
    print(f"  Expected one of: {sp2_electrophiles}")
    print(f"  But only detected: {list(reactant_profile.keys())}")

print("\nChecking for organoboron:")
found_organoboron = [m for m in organoboron if m in reactant_profile]
if found_organoboron:
    print(f"  ✓ FOUND: {found_organoboron}")
else:
    print("  ✗ NOT FOUND")
    print(f"  Expected one of: {organoboron}")

print("\nChecking for coupling products:")
found_products = [m for m in coupling_products if m in product_profile]
if found_products:
    print(f"  ✓ FOUND: {found_products}")
else:
    print("  ✗ NOT FOUND")
    print(f"  Expected one of: {coupling_products}")
    print(f"  But only detected: {list(product_profile.keys())}")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)

if found_electrophiles and found_organoboron:
    print("✓ Reactant motifs match Suzuki requirements!")
    if not found_products:
        print("✗ However, product motifs don't match (but this might not be checked)")
else:
    print("✗ Suzuki detection failing due to missing reactant motifs:")
    if not found_electrophiles:
        print("  - Missing sp2_electrophile")
    if not found_organoboron:
        print("  - Missing organoboron nucleophile")
