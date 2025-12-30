
from chemtools.featurizers.motif_registry import build_compound_detect_registry, _default_registry_paths
from chemtools.util.rdkit_helpers import parse_smiles
from rdkit import Chem

paths = _default_registry_paths()
registry = build_compound_detect_registry(paths)

compiled = registry.get("compiled_compounds", {})
if "Any-OH" in compiled:
    print(f"Any-OH queries: {len(compiled['Any-OH'])}")
    # We can't easily print the query object, but we can test it
    mol = parse_smiles("OCC1CCCCC1")
    for query in compiled["Any-OH"]:
        match = mol.HasSubstructMatch(query)
        print(f"Match OCC1CCCCC1: {match}")
else:
    print("Any-OH not found in registry")

# Check what the SMARTS was
import json
from pathlib import Path
with open(paths["compounds"], "r") as f:
    compounds = json.load(f)["compounds"]
    any_oh = next(c for c in compounds if c["id"] == "Any-OH")
    print(f"Any-OH entry: {any_oh}")

with open(paths["groups"], "r") as f:
    groups = json.load(f)["groups"]
    r_group = next(g for g in groups if g["id"] == "R")
    oh_group = next(g for g in groups if g["id"] == "OH")
    print(f"R group: {r_group}")
    print(f"OH group: {oh_group}")
