import traceback
from chemtools import detect_reaction_type

try:
    r = detect_reaction_type('CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1')
    print(f'Match: {r.top_match}')
    print(f'Error: {r.error}')
    if r.top_match:
        print(f'{r.top_match.reaction_type}: {r.top_match.electrophile} + {r.top_match.nucleophile} -> {r.top_match.product}')
except Exception as e:
    traceback.print_exc()
