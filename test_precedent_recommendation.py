"""
Test precedent-based recommendation for a specific reaction.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

from chemtools import precedent
from chemtools.featurizers import reaction_pair as feat_pair
from chemtools.recommend.utils import pick_electrophile_nucleophile
from chemtools.smiles import normalize_reaction

# Test reaction SMILES
reaction_smiles = "c1ccc2c(c1)CCNC2.Nc1ccc(F)c(F)c1>>Nc1ccc(F)c(N2CCc3ccccc3C2)c1"

print(f"Testing precedent recommendation for:")
print(f"Reaction: {reaction_smiles}\n")

# Parse reactants
normalized = normalize_reaction(reaction_smiles)
reactant_pool = []
for entry in normalized.get("reactants", []) or []:
    if not isinstance(entry, dict):
        continue
    smi = entry.get("smiles_norm") or entry.get("largest_smiles") or entry.get("input")
    if smi:
        reactant_pool.append(smi)

print(f"Reactants: {reactant_pool}")

# Pick electrophile and nucleophile
elec, nuc = pick_electrophile_nucleophile(reactant_pool)
print(f"Electrophile: {elec}")
print(f"Nucleophile: {nuc}\n")

# Featurize the pair
features = feat_pair.featurize_pair(elec, nuc).get("flat", {}) if (elec or nuc) else {}
print(f"Features extracted: {len(features)} features\n")

# Call precedent search with DRFP
relax = {
    "use_drfp": True,
    "reaction_smiles": reaction_smiles,
    "filter_by_reagent_database": False,
}

print("Running precedent search (top 10)...")
print("=" * 80)

pack = precedent.knn(family=None, features=features, k=10, relax=relax)
precedents = list(pack.get("precedents", []) or [])

print(f"\nFound {len(precedents)} precedent recommendations:\n")

for i, prec in enumerate(precedents, 1):
    conditions = prec.get("conditions") or {}
    
    catalyst = conditions.get("catalyst", "")
    ligand = conditions.get("ligand", "")
    base = conditions.get("base", prec.get("base_uid", ""))
    solvent = conditions.get("solvent", "")
    additive = conditions.get("additive", "")
    
    yield_val = prec.get("yield") or conditions.get("yield_pct")
    temp = conditions.get("temperature_c")
    time_h = conditions.get("time_h")
    
    rxn_type = prec.get("rxn_type", "")
    rxn_id = prec.get("reaction_id", "")
    similarity = prec.get("similarity", 0.0)
    
    print(f"Rank {i}:")
    print(f"  Similarity: {similarity:.3f}")
    print(f"  Reaction Type: {rxn_type}")
    print(f"  Reaction ID: {rxn_id}")
    print(f"  Catalyst: {catalyst or 'None'}")
    print(f"  Ligand: {ligand or 'None'}")
    print(f"  Base: {base or 'None'}")
    print(f"  Solvent: {solvent or 'None'}")
    print(f"  Additive: {additive or 'None'}")
    if yield_val is not None:
        print(f"  Yield: {yield_val}%")
    if temp is not None:
        print(f"  Temperature: {temp}°C")
    if time_h is not None:
        print(f"  Time: {time_h}h")
    print()

print("=" * 80)
print("Done!")
