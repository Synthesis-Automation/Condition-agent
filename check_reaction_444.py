"""
Check what's in the precedent database for reaction IDs 443 and 444.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from chemtools import precedent

# Load all precedents
pack = precedent.knn(family=None, features={}, k=100000, relax={})
precedents = list(pack.get("precedents", []) or [])

print(f"Total precedents loaded: {len(precedents)}\n")

# Look for reactions 443 and 444
target_ids = ["C-N Coupling:443", "C-N Coupling:444", "C-N Coupling:442"]

for target_id in target_ids:
    found = [p for p in precedents if p.get("reaction_id") == target_id]
    if found:
        prec = found[0]
        print(f"Found: {target_id}")
        print(f"  Reaction SMILES: {prec.get('reaction_smiles', 'N/A')}")
        conditions = prec.get("conditions") or {}
        print(f"  Catalyst: {conditions.get('catalyst', 'None')}")
        print(f"  Ligand: {conditions.get('ligand', 'None')}")
        print(f"  Base: {conditions.get('base', prec.get('base_uid', 'None'))}")
        print(f"  Solvent: {conditions.get('solvent', prec.get('solvent_uid', 'None'))}")
        print(f"  Additive: {conditions.get('additive', 'None')}")
        print(f"  Yield: {prec.get('yield', 'N/A')}%")
        print(f"  Reaction Type: {prec.get('rxn_type', 'N/A')}")
        print()
    else:
        print(f"NOT FOUND: {target_id}\n")
