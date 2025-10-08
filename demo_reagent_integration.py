"""
Demo: ChemTools v2.0 with Reagent Lookup Integration

Shows the complete ChemTools v2.0 API including the new reagent namespace.
"""

from chemtools import chem

print("=" * 70)
print("ChemTools v2.0 - Reagent Lookup Demo")
print("=" * 70)

# ============================================================================
# 1. List Available Reagent Databases
# ============================================================================
print("\n📚 Available Reagent Databases:")
print("-" * 70)
types = chem.reagent.list_types()
for reagent_type in types:
    print(f"  • {reagent_type}")

# ============================================================================
# 2. Find Specific Reagents
# ============================================================================
print("\n🔍 Finding Reagents:")
print("-" * 70)

# Find ligand by abbreviation
ligand = chem.reagent.find("PPh3", "ligand")
if ligand:
    print(f"✓ Found: {ligand['name']}")
    print(f"  CAS: {ligand.get('cas', 'N/A')}")
    print(f"  Abbreviations: {ligand.get('abbreviation', [])}")

# Find base by CAS number
base = chem.reagent.find("584-08-7", "base")
if base:
    print(f"\n✓ Found: {base['name']}")
    print(f"  Abbreviation: {base.get('abbreviation', ['N/A'])[0]}")
    print(f"  SMILES: {base.get('smiles', 'N/A')}")

# Find solvent by name
solvent = chem.reagent.find("dimethylformamide", "solvent")
if solvent:
    print(f"\n✓ Found: {solvent['name']}")
    print(f"  Aliases: {solvent.get('aliases', [])[:3]}")

# ============================================================================
# 3. Enrich Reagent Information
# ============================================================================
print("\n\n📝 Enriching Reagent Information:")
print("-" * 70)

# Enrich a metal precursor
catalyst_info = chem.reagent.enrich("Pd(PPh3)4", "metal_precursor")
print(f"Catalyst: {catalyst_info.get('name', 'Not found')}")
print(f"Found in DB: {catalyst_info.get('found', False)}")
if catalyst_info.get('found'):
    print(f"CAS: {catalyst_info.get('cas', 'N/A')}")
    print(f"SMILES: {catalyst_info.get('smiles', 'N/A')[:50]}...")

# ============================================================================
# 4. Enrich Complete Condition Set
# ============================================================================
print("\n\n🧪 Enriching Full Condition Set:")
print("-" * 70)

# Example conditions from a recommendation
conditions = {
    "catalyst": "Pd(OAc)2",
    "ligand": "XPhos",
    "base": "K2CO3",
    "solvent": "DMF",
    "temperature": "100 °C"
}

print("Original conditions:")
for key, value in conditions.items():
    print(f"  {key}: {value}")

# Enrich with reagent details
enriched = chem.reagent.enrich_conditions(conditions)

print("\nEnriched conditions (with _details fields):")
for key in conditions.keys():
    details_key = f"{key}_details"
    if details_key in enriched:
        details = enriched[details_key]
        if details.get('found'):
            print(f"  {key}_details:")
            print(f"    ✓ Name: {details.get('name', 'N/A')}")
            print(f"    ✓ CAS: {details.get('cas', 'N/A')}")
        else:
            print(f"  {key}_details: Not found in database")

# ============================================================================
# 5. Integration with Other ChemTools Features
# ============================================================================
print("\n\n🔗 Full ChemTools v2.0 API:")
print("-" * 70)

print("Core operations:")
print("  ✓ chem.smiles.normalize()          - SMILES normalization")
print("  ✓ chem.router.detect_family()      - Reaction family detection")
print("  ✓ chem.properties.lookup()         - Compound properties")
print("  ✓ chem.constraints.filter()        - Constraint filtering")
print("  ✓ chem.reagent.find()              - Reagent lookup (NEW!)")

print("\nData operations:")
print("  ✓ chem.precedent.knn()             - Precedent search")
print("  ✓ chem.recommend.conditions()      - ML recommendations")
print("  ✓ chem.explain.precedents()        - Generate explanations")

print("\n" + "=" * 70)
print("✅ ChemTools v2.0 with Reagent Lookup Ready!")
print("=" * 70)
