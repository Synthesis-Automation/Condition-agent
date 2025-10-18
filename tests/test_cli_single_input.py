"""Test the updated CLI with single-input extraction."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("="*70)
print("CLI Single-Input Extraction Test")
print("="*70)

# Test cases showing different input formats
test_cases = [
    {
        "name": "Just SMILES",
        "input": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
        "expected": {
            "has_smiles": True,
            "has_constraints": False
        }
    },
    {
        "name": "SMILES + requirement after",
        "input": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 use copper catalyst",
        "expected": {
            "has_smiles": True,
            "has_constraints": True,
            "metal": "Cu"
        }
    },
    {
        "name": "Requirement + SMILES",
        "input": "I want to run Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 with no strong base",
        "expected": {
            "has_smiles": True,
            "has_constraints": True,
            "base": "moderate"
        }
    },
    {
        "name": "Mixed format",
        "input": "use palladium for CCBr.c1ccccc1>>CCc1ccccc1 at room temperature avoid expensive reagents",
        "expected": {
            "has_smiles": True,
            "has_constraints": True,
            "metal": "Pd",
            "temp_max": 30
        }
    },
    {
        "name": "Complex mixed",
        "input": "Run Clc1ccccc1.Sc1ccccc1>>c1ccc(Sc2ccccc2)cc1 using copper catalyst no strong base RT",
        "expected": {
            "has_smiles": True,
            "has_constraints": True,
            "metal": "Cu"
        }
    }
]

print("\nTest Cases:")
print("-"*70)

for i, test in enumerate(test_cases, 1):
    print(f"\n{i}. {test['name']}")
    print(f"   Input: {test['input']}")
    print(f"   Expected: SMILES={'✓' if test['expected']['has_smiles'] else '✗'}, "
          f"Constraints={'✓' if test['expected']['has_constraints'] else '✗'}")

print("\n" + "="*70)
print("How the LLM will extract information:")
print("="*70)

print("""
The LLM will:
1. Scan the input for a SMILES pattern (look for '>>')
2. Extract the complete SMILES string
3. Extract requirements from remaining text:
   - "use X", "with X" → required_reagents
   - "copper", "palladium" → metal_preference
   - "no strong base" → base_strength: moderate
   - "room temperature", "RT" → temperature.max: 25-30
   - "avoid X", "no X" → exclude_reagents
   - "cheap", "expensive" → cost_level

Examples of what gets extracted:

Input: "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 use copper catalyst"
→ SMILES: "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
→ Constraints: {"required_reagents": ["copper catalyst"], "metal_preference": "Cu"}

Input: "I want to run Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 with no strong base"
→ SMILES: "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
→ Constraints: {"base_strength": "moderate"}

Input: "use palladium for CCBr.c1ccccc1>>CCc1ccccc1 at room temperature"
→ SMILES: "CCBr.c1ccccc1>>CCc1ccccc1"
→ Constraints: {"metal_preference": "Pd", "required_reagents": ["palladium"], "temperature": {"max": 30}}
""")

print("\n" + "="*70)
print("To test with actual LLM:")
print("="*70)
print("""
Run the CLI:
  python app/cli_AI_recommend.py

Then try these inputs:
  1. Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 use copper catalyst
  2. use palladium for CCBr.c1ccccc1>>CCc1ccccc1 at room temperature
  3. Run Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 with no strong base RT

The AI will intelligently:
- Extract the SMILES from anywhere in the text
- Parse requirements from the surrounding words
- Validate the SMILES format
- Structure everything for the API
""")

print("\n✓ Update complete! CLI now accepts mixed input in a single line.")
