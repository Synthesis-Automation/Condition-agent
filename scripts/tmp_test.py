import os, sys
sys.path.insert(0, os.path.abspath('.'))
from chemtools import smiles
s='c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1'
print('normalize_reaction output:')
res = smiles.normalize_reaction(s)
for k in ('input','normalized','errors'):
    print(k, ':', res.get(k))
print('reactants:', [x.get('smiles_norm') or x.get('largest_smiles') or x.get('input') for x in res.get('reactants')])
print('products:', [x.get('smiles_norm') or x.get('largest_smiles') or x.get('input') for x in res.get('products')])
