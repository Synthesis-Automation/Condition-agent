#!/usr/bin/env python
"""Debug: Why does coupling confirmation fail for C_O_Coupling?"""

from chemtools.featurizers.detection.coupling import (
    get_coupling_confirmation_specs,
    _compile_mapped_patterns,
    _find_attachment_sites,
)
from chemtools.util.rdkit_helpers import parse_smiles

# Phenol + iodobenzonitrile → aryl ether
rxn = 'Cc1cc(C)cc(O)c1.N#Cc1ccc(I)cc1>>Cc1cc(C)cc(Oc2ccc(C#N)cc2)c1'
reactants_str, products_str = rxn.split('>>')
reactants = [parse_smiles(s) for s in reactants_str.split('.')]

print("Testing coupling confirmation for C_O_Coupling...")
print(f"Reactant 1: {reactants_str.split('.')[0]} (phenol)")
print(f"Reactant 2: {reactants_str.split('.')[1]} (iodobenzonitrile)")
print()

# Get specs
specs = get_coupling_confirmation_specs()
spec = specs.get('C_O_Coupling')

print(f"Electrophile SMARTS: {spec.electrophile_smarts}")
print(f"Nucleophile SMARTS: {spec.nucleophile_smarts}")
print()

# Find sites
electrophile_sites = []
nucleophile_sites = []
for idx, mol in enumerate(reactants):
    print(f"Checking reactant {idx}...")
    e_sites = _find_attachment_sites(mol, spec.electrophile_smarts)
    n_sites = _find_attachment_sites(mol, spec.nucleophile_smarts)
    print(f"  Electrophile sites: {e_sites}")
    print(f"  Nucleophile sites: {n_sites}")
    
    for attach_idx, lg_idx in e_sites:
        if lg_idx is not None:
            electrophile_sites.append((idx, attach_idx, lg_idx))
    for attach_idx, lg_idx in n_sites:
        nucleophile_sites.append((idx, attach_idx, lg_idx))

print()
print(f"Total electrophile sites: {len(electrophile_sites)}")
print(f"Total nucleophile sites: {len(nucleophile_sites)}")
print()

if len(electrophile_sites) != 1:
    print(f"❌ Problem: Expected 1 electrophile site, found {len(electrophile_sites)}")
if len(nucleophile_sites) != 1:
    print(f"❌ Problem: Expected 1 nucleophile site, found {len(nucleophile_sites)}")
