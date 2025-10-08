"""Quick test of rxn-insight detector."""

import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from chemtools.reaction_type_detector import detect_reaction_type

smiles = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"

result = detect_reaction_type(smiles)

print(f"Available: {result['available']}")
print(f"Success: {result['success']}")
print(f"Class: {result['rxn_class']}")
print(f"Name: {result['rxn_name']}")
print(f"Mapped: {result['mapped_family']}")
print(f"Catalysts: {result['catalysts']}")
print(f"Confidence: {result['confidence']}")
