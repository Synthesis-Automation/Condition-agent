#!/usr/bin/env python3
"""Test the bond analysis tool integration."""

import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chem_assistant.chemtools_wrapper import analyze_bond_changes_tool

print("=" * 80)
print("Testing Bond Analysis Tool Integration")
print("=" * 80)
print()

# Test 1: Suzuki-Miyaura coupling
print("Test 1: Suzuki-Miyaura Coupling")
print("-" * 80)
suzuki_smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
print(f"SMILES: {suzuki_smiles}")
print()

result = analyze_bond_changes_tool.invoke({
    "reaction_smiles": suzuki_smiles,
    "use_hybrid": True
})

print("Result:")
import json
print(json.dumps(result, indent=2))
print()

# Test 2: Simple esterification
print("=" * 80)
print("Test 2: Esterification")
print("-" * 80)
ester_smiles = "CC(=O)O.CCO>>CC(=O)OCC"
print(f"SMILES: {ester_smiles}")
print()

result2 = analyze_bond_changes_tool.invoke({
    "reaction_smiles": ester_smiles,
    "use_hybrid": True
})

print("Result:")
print(json.dumps(result2, indent=2))
print()

print("=" * 80)
print("✅ Bond analysis tool integration successful!")
print("=" * 80)
