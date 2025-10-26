"""Quick test to show reactant type information"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chemtools.analysis import analyze_reaction
import json

# Test a Buchwald-Hartwig reaction
smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
print("Testing reaction:", smiles)
print("Description: Bromobenzene + Aniline -> Diphenylamine (Buchwald-Hartwig)")
print("\n" + "="*80)

result = analyze_reaction(smiles)

print(f"\nReaction Family: {result['family']['canonical_id']}")
print(f"Number of Reactants: {len(result['reactants'])}")

# Print full structure to see what's available
print("\n" + "="*80)
print("Full Analysis Structure:")
print(json.dumps(result, indent=2, default=str))
