#!/usr/bin/env python3
"""
Demo: Bond Breaking and Formation Analysis with Atom-Mapped SMILES
===================================================================

This demonstrates the EXISTING functionality for analyzing which bonds
break and form in a reaction, using atom-mapped SMILES.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

# Check if RDKit is available
try:
    from chemtools.util.reaction_center_detector import identify_changed_atoms_from_mapped_smiles
    print("✅ reaction_center_detector imported successfully")
except ImportError as e:
    print(f"❌ Could not import reaction_center_detector: {e}")
    sys.exit(1)

print("="*80)
print("BOND ANALYSIS DEMO - Atom-Mapped Reactions")
print("="*80)

# Example 1: Sonogashira Coupling
print("\n📌 Example 1: Sonogashira Coupling")
print("-" * 80)

# Manually atom-mapped Sonogashira: Ar-I + HC≡C-R → Ar-C≡C-R
sonogashira_mapped = "[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[I:7].[C:8]#[C:9][CH3:10]>>[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[C:8]#[C:9][CH3:10]"

print(f"Reaction: {sonogashira_mapped}")
print()

result1 = identify_changed_atoms_from_mapped_smiles(sonogashira_mapped)

if result1.get('success'):
    print("✅ Analysis successful!")
    print()
    print(f"Changed atoms (map numbers): {sorted(result1['changed_atoms'])}")
    print(f"  → Atoms 6, 7, 8 are involved in bond changes")
    print()
    print(f"Broken bonds: {result1['broken_bonds']}")
    print(f"  → Bond between atoms 6-7 (Ar-I bond) is BROKEN")
    print()
    print(f"Formed bonds: {result1['formed_bonds']}")
    print(f"  → Bond between atoms 6-8 (Ar-C bond) is FORMED")
    print()
    print(f"Spectator atoms: {sorted(result1['spectator_atoms'])}")
    print(f"  → Atoms {sorted(result1['spectator_atoms'])} don't change connectivity")
else:
    print(f"❌ Error: {result1.get('error')}")

# Example 2: Suzuki-Miyaura Coupling
print("\n" + "="*80)
print("📌 Example 2: Suzuki-Miyaura Coupling")
print("-" * 80)

# Ar-Br + Ar'-B(OH)2 → Ar-Ar'
suzuki_mapped = "[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[Br:7].[c:8]1[cH:9][cH:10][cH:11][cH:12][c:13]1[B:14]([OH:15])[OH:16]>>[c:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[c:8]1[cH:9][cH:10][cH:11][cH:12][c:13]1"

print(f"Reaction: Ar-Br + Ar'-B(OH)2 → Ar-Ar'")
print()

result2 = identify_changed_atoms_from_mapped_smiles(suzuki_mapped)

if result2.get('success'):
    print("✅ Analysis successful!")
    print()
    print(f"Changed atoms: {sorted(result2['changed_atoms'])}")
    print()
    print(f"Broken bonds: {result2['broken_bonds']}")
    print(f"  → Ar-Br bond broken (6-7)")
    print(f"  → Ar'-B bonds broken (13-14, 14-15, 14-16)")
    print()
    print(f"Formed bonds: {result2['formed_bonds']}")
    print(f"  → New Ar-Ar' bond formed (6-13)")
    print()
    spectators = result2['spectator_atoms']
    print(f"Spectator atoms ({len(spectators)}): Most aromatic carbons unchanged")
else:
    print(f"❌ Error: {result2.get('error')}")

# Example 3: Simple Amide Coupling
print("\n" + "="*80)
print("📌 Example 3: Amide Coupling")
print("-" * 80)

# R-COOH + R'-NH2 → R-CO-NH-R' (simplified, condensation not shown)
amide_mapped = "[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>>[CH3:1][C:2](=[O:3])[NH:5][CH3:6]"

print(f"Reaction: R-COOH + R'-NH2 → R-CO-NH-R'")
print()

result3 = identify_changed_atoms_from_mapped_smiles(amide_mapped)

if result3.get('success'):
    print("✅ Analysis successful!")
    print()
    print(f"Changed atoms: {sorted(result3['changed_atoms'])}")
    print()
    print(f"Broken bonds: {result3['broken_bonds']}")
    print(f"  → C-OH bond broken (2-4)")
    print()
    print(f"Formed bonds: {result3['formed_bonds']}")
    print(f"  → C-N amide bond formed (2-5)")
    print()
else:
    print(f"❌ Error: {result3.get('error')}")

# Summary
print("\n" + "="*80)
print("📊 SUMMARY")
print("="*80)
print("""
✅ What works:
   • Bond breaking/formation detection with atom-mapped SMILES
   • Identifies exact atoms involved in changes
   • Distinguishes reaction center from spectator atoms

❌ What's missing:
   • Automatic atom mapping for unmapped SMILES
   • RXNMapper integration (would enable this)
   • Public API exposure (currently in util/, not main API)

💡 To use with unmapped reactions:
   1. Install RXNMapper: pip install rxnmapper
   2. Add mapping: from rxnmapper import RXNMapper; mapper.get_attention_guided_atom_maps([smiles])
   3. Then use reaction_center_detector on mapped result
""")

print("="*80)
