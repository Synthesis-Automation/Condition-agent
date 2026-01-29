import sys
# Force reimport
if 'chemtools' in sys.modules:
    to_reload = [k for k in sys.modules.keys() if k.startswith('chemtools')]
    for mod in to_reload:
        del sys.modules[mod]

from chemtools.featurizers.formatters.reaction import featurize_reaction

rxn = 'CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1'

print("=" * 70)
print("Dithiocarbamate Motif Detection Test (after adding Thiocarbonyl-SH)")
print("=" * 70)

result = featurize_reaction(rxn)

print("\nREACTANT MOTIFS:")
for i, r in enumerate(result['reactants']):
    motifs = sorted(set(m.get('id') or m.get('compound_id', 'unknown') for m in r.get('motifs', [])))
    print(f"  Reactant {i}: {motifs}")

thiol_found = False
for r in result['reactants']:
    for m in r.get('motifs', []):
        mid = m.get('id') or m.get('compound_id', '')
        if 'SH' in mid:
            thiol_found = True
            print(f"\n✓ THIOL FOUND: {mid}")
            break
    if thiol_found:
        break

if not thiol_found:
    print("\n❌ NO THIOL MOTIF DETECTED")
    print("\nReactant 0 (dithiocarbamate) details:")
    if result['reactants']:
        r0 = result['reactants'][0]
        print(f"  SMILES: {r0.get('smiles')}")
        print(f"  Motifs: {[m.get('id') or m.get('compound_id') for m in r0.get('motifs', [])]}")
