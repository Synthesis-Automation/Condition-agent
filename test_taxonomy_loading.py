#!/usr/bin/env python3
"""Test if taxonomy is loading correctly."""

from chemtools.taxonomy import load_registry

registry = load_registry()

print(f"Reactant types loaded: {len(list(registry.iter_reactant_types()))}")
print(f"Reaction types loaded: {len(list(registry.iter_reaction_types()))}")

# Test a few specific types
arx = registry.get_reactant_type("arx")
print(f"\nArX reactant type: {arx.name if arx else 'NOT FOUND'}")

epoxide = registry.get_reactant_type("epoxide")
print(f"Epoxide reactant type: {epoxide.name if epoxide else 'NOT FOUND'}")

diene = registry.get_reactant_type("diene")
print(f"Diene reactant type: {diene.name if diene else 'NOT FOUND'}")

# Test reaction types
diels_alder = registry.get_reaction_type("diels_alder")
print(f"\nDiels-Alder reaction: {diels_alder.name if diels_alder else 'NOT FOUND'}")

hydrogenation = registry.get_reaction_type("hydrogenation")
print(f"Hydrogenation reaction: {hydrogenation.name if hydrogenation else 'NOT FOUND'}")
