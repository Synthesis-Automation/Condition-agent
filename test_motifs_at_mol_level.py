import sys
# Force reimport
if 'chemtools' in sys.modules:
    to_reload = [k for k in sys.modules.keys() if k.startswith('chemtools')]
    for mod in to_reload:
        del sys.modules[mod]

from chemtools.context import ReactionContext

rxn = "CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1"

print("=== Checking Motif Detection at Molecular Level ===\n")

ctx = ReactionContext.from_smiles(rxn)

print("Reactant 0 (dithiocarbamate):")
if ctx.reactants and len(ctx.reactants) > 0:
    r0_motifs = ctx.reactants[0].motifs
    print(f"  Motifs ({len(r0_motifs)}):")
    for m in r0_motifs:
        mid = m.get('compound_id') or m.get('id', 'unknown')
        print(f"    - {mid}")
    
    # Check for thiol
    has_thiol = any('SH' in (m.get('compound_id') or m.get('id', '')) for m in r0_motifs)
    print(f"  Has thiol motif: {has_thiol}")
    
print("\nReactant 2 (aryl iodide):")
if ctx.reactants and len(ctx.reactants) > 2:
    r2_motifs = ctx.reactants[2].motifs
    print(f"  Motifs ({len(r2_motifs)}):")
    for m in r2_motifs:
        mid = m.get('compound_id') or m.get('id', 'unknown')
        print(f"    - {mid}")
        
print("\nProduct (aryl sulfide):")
if ctx.products and len(ctx.products) > 0:
    p0_motifs = ctx.products[0].motifs
    print(f"  Motifs ({len(p0_motifs)}):")
    for m in p0_motifs:
        mid = m.get('compound_id') or m.get('id', 'unknown')
        print(f"    - {mid}")
