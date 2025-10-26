#!/usr/bin/env python
"""Debug why Heck/Kumada/Negishi aren't being detected."""

import sys
sys.path.insert(0, '.')

from chemtools.router import detect_family, _rule_hits

# Test Heck reaction
heck_reactants = ["Brc1ccccc1", "C=Cc1ccccc1"]
print("=" * 80)
print("Testing Heck detection directly")
print("=" * 80)
print(f"Reactants: {heck_reactants}")

# Check hits
hits = _rule_hits(heck_reactants)
print(f"\nHits detected:")
for k, v in hits.items():
    if v:
        print(f"  {k}: {v}")

# Check detect_family
result = detect_family(heck_reactants)
print(f"\ndetect_family result:")
print(f"  family: {result.get('family')}")
print(f"  confidence: {result.get('confidence')}")
print(f"  hits: {result.get('hits')}")

# Test Kumada
print("\n" + "=" * 80)
print("Testing Kumada detection directly")
print("=" * 80)
kumada_reactants = ["Brc1ccccc1", "C[Mg]Br"]
print(f"Reactants: {kumada_reactants}")

hits = _rule_hits(kumada_reactants)
print(f"\nHits detected:")
for k, v in hits.items():
    if v:
        print(f"  {k}: {v}")

result = detect_family(kumada_reactants)
print(f"\ndetect_family result:")
print(f"  family: {result.get('family')}")
print(f"  confidence: {result.get('confidence')}")

# Test Negishi
print("\n" + "=" * 80)
print("Testing Negishi detection directly")
print("=" * 80)
negishi_reactants = ["Brc1ccccc1", "C[Zn]Cl"]
print(f"Reactants: {negishi_reactants}")

hits = _rule_hits(negishi_reactants)
print(f"\nHits detected:")
for k, v in hits.items():
    if v:
        print(f"  {k}: {v}")

result = detect_family(negishi_reactants)
print(f"\ndetect_family result:")
print(f"  family: {result.get('family')}")
print(f"  confidence: {result.get('confidence')}")
