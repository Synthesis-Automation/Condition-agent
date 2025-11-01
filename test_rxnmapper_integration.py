#!/usr/bin/env python3
"""
Demo: RXNMapper Integration for Bond Analysis
==============================================

This demonstrates the NEW RXNMapper integration that enables:
1. Automatic atom mapping for unmapped reactions
2. Bond breaking/formation analysis
3. Reaction center identification
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

print("="*80)
print("RXNMAPPER INTEGRATION DEMO")
print("="*80)

# Check availability
print("\n1️⃣  Checking RXNMapper availability...")
print("-" * 80)

from chemtools import rxnmapper_available

if rxnmapper_available():
    print("✅ RXNMapper is AVAILABLE")
else:
    print("❌ RXNMapper is NOT available")
    print("   Install with: pip install rxnmapper")
    print("\n⚠️  Continuing with demo (will show error messages)")

# Test automatic atom mapping
print("\n2️⃣  Testing Automatic Atom Mapping...")
print("-" * 80)

from chemtools import add_atom_mapping

# Unmapped Suzuki-Miyaura coupling
unmapped_suzuki = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

print(f"Input (unmapped): {unmapped_suzuki}")
print()

mapping_result = add_atom_mapping(unmapped_suzuki)

if mapping_result['success']:
    print("✅ Atom mapping succeeded!")
    print(f"Method: {mapping_result['method']}")
    if mapping_result.get('confidence'):
        print(f"Confidence: {mapping_result['confidence']:.3f}")
    print(f"\nMapped SMILES:")
    print(f"  {mapping_result['mapped_smiles'][:100]}...")
else:
    print(f"❌ Atom mapping failed: {mapping_result.get('error')}")

# Test high-level bond analysis
print("\n3️⃣  Testing High-Level Bond Analysis...")
print("-" * 80)

from chemtools import analyze_bond_changes

# Unmapped Sonogashira coupling
unmapped_sono = "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"

print(f"Reaction: Sonogashira coupling")
print(f"Input: {unmapped_sono}")
print()

bond_result = analyze_bond_changes(unmapped_sono, auto_map=True)

if bond_result['success']:
    print("✅ Bond analysis succeeded!")
    print()
    print(f"Changed atoms: {sorted(bond_result['changed_atoms'])}")
    print(f"Broken bonds: {bond_result['broken_bonds']}")
    print(f"  → Ar-Br bond breaks")
    print(f"Formed bonds: {bond_result['formed_bonds']}")
    print(f"  → Ar-C≡C bond forms")
    print(f"Spectator atoms: {len(bond_result['spectator_atoms'])} atoms unchanged")
    
    if bond_result.get('mapping_confidence'):
        print(f"\nMapping confidence: {bond_result['mapping_confidence']:.3f}")
else:
    print(f"❌ Bond analysis failed: {bond_result.get('error')}")

# Test with multiple reactions
print("\n4️⃣  Testing Batch Processing...")
print("-" * 80)

from chemtools._atom_mapping import batch_add_atom_mapping

reactions = [
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",  # Buchwald-Hartwig
    "Ic1ccccc1.C#C>>C(#C)c1ccccc1",  # Sonogashira
    "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",  # Suzuki
]

print(f"Processing {len(reactions)} reactions in batch...")
print()

batch_results = batch_add_atom_mapping(reactions)

for i, (rxn, result) in enumerate(zip(reactions, batch_results), 1):
    if result['success']:
        print(f"  {i}. ✅ Mapped successfully")
        if result.get('confidence'):
            print(f"     Confidence: {result['confidence']:.3f}")
    else:
        print(f"  {i}. ❌ Failed: {result.get('error')}")

# Test with compare_unmapped_reaction_to_find_changes
print("\n5️⃣  Testing Auto-Map + Analyze (Single Function)...")
print("-" * 80)

from chemtools import compare_unmapped_reaction_to_find_changes

# Unmapped amide coupling
unmapped_amide = "CC(=O)O.CCN>>CC(=O)NCC"

print(f"Reaction: Simple amide coupling")
print(f"Input: {unmapped_amide}")
print()

auto_result = compare_unmapped_reaction_to_find_changes(unmapped_amide)

if auto_result['success']:
    print("✅ Auto-mapping and analysis succeeded!")
    print()
    print(f"Broken bonds: {auto_result['broken_bonds']}")
    print(f"  → C-OH bond breaks")
    print(f"Formed bonds: {auto_result['formed_bonds']}")
    print(f"  → C-N amide bond forms")
    
    if auto_result.get('auto_mapped'):
        print(f"\n✨ Automatically mapped with confidence: {auto_result.get('mapping_confidence', 'N/A')}")
else:
    print(f"❌ Analysis failed: {auto_result.get('error')}")
    if auto_result.get('message'):
        print(f"   {auto_result['message']}")
    if auto_result.get('recommendation'):
        print(f"   💡 {auto_result['recommendation']}")

# Integration with detect_reaction
print("\n6️⃣  Integration with Reaction Detection...")
print("-" * 80)

from chemtools import detect_reaction

print("Detecting reaction type AND analyzing bonds...")
print()

detection = detect_reaction(unmapped_suzuki, use_ml=True)
bonds = analyze_bond_changes(unmapped_suzuki, auto_map=True)

print(f"Detected reaction type: {detection['family']}")
print(f"Detection confidence: {detection['confidence']:.2f}")
print()

if bonds['success']:
    print(f"Bond analysis:")
    print(f"  Broken: {len(bonds['broken_bonds'])} bonds")
    print(f"  Formed: {len(bonds['formed_bonds'])} bonds")
    print(f"  Changed atoms: {len(bonds['changed_atoms'])}")

# Summary
print("\n" + "="*80)
print("📊 SUMMARY: NEW CAPABILITIES")
print("="*80)
print("""
✅ What's NEW:

1. Automatic Atom Mapping
   • add_atom_mapping(reaction_smiles) → mapped SMILES
   • batch_add_atom_mapping([reactions]) → batch mapping
   • Uses RXNMapper attention-guided mapping

2. High-Level Bond Analysis
   • analyze_bond_changes(reaction_smiles, auto_map=True)
   • Automatically maps if needed, then analyzes
   • Returns broken_bonds, formed_bonds, changed_atoms

3. Unified API
   • compare_unmapped_reaction_to_find_changes(smiles)
   • One-step: unmapped SMILES → bond analysis
   • Graceful degradation if RXNMapper not installed

4. Public API Exposure
   • All functions available via: from chemtools import ...
   • Documented in chemtools.__all__
   • Type hints and docstrings

🔧 Usage Patterns:

   # Simple: Just get mapped SMILES
   result = add_atom_mapping("Br.C#C>>C#C")
   
   # Advanced: Analyze bonds (auto-maps if needed)
   bonds = analyze_bond_changes("Br.C#C>>C#C", auto_map=True)
   
   # Combined: Detection + Bond Analysis
   rxn_type = detect_reaction(smiles)
   bonds = analyze_bond_changes(smiles)

📦 Installation:

   pip install rxnmapper

   Note: RXNMapper is optional. If not installed, functions will
   gracefully return error messages with installation instructions.
""")

print("="*80)
