
from chemtools.taxonomy.reaction_catalog import load_reaction_catalog

defs, aliases = load_reaction_catalog()

target_ids = [
    "oxidation_alcohol_to_carbonyl",
    "reduction_carbonyl_to_alcohol",
    "reduction_nitro_to_amine",
    "halogenation_aromatic"
]

for tid in target_ids:
    if tid in defs:
        d = defs[tid]
        print(f"ID: {tid}")
        print(f"  Reactants: {d.reactants}")
    else:
        print(f"ID: {tid} NOT FOUND")
