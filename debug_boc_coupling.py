#!/usr/bin/env python
"""Debug Boc-amine + aryl bromide C-N coupling detection."""

from chemtools.featurizers.formatters.reaction import featurize_reaction
from chemtools.reaction_type_detection import detect_reaction_types

rxn = 'O=c1oc2cc(Br)ccc2cc1-c1ccccc1.CC(C)(C)OC(=O)NCCN>>NCCNc1ccc2cc(-c3ccccc3)c(=O)oc2c1'

print("=" * 70)
print("Reaction: Boc-protected amine + aryl bromide → aromatic amine")
print("=" * 70)

result = featurize_reaction(rxn)

print("\n1. REACTANT MOTIFS:")
for i, r in enumerate(result['reactants']):
    motifs = sorted(set(m['id'] for m in r.get('motifs', [])))
    print(f"   Reactant {i}: {motifs}")

print("\n2. PRODUCT MOTIFS:")
for i, p in enumerate(result['products']):
    motifs = sorted(set(m['id'] for m in p.get('motifs', [])))
    print(f"   Product {i}: {motifs}")

print(f"\n3. REACTION KEY:")
print(f"   {result.get('reaction_key', 'N/A')}")

agg = result.get('aggregates', {})
print(f"\n4. MOTIF CLASSIFICATION:")
print(f"   Reacted: {sorted(agg.get('reacted', []))}")
print(f"   Formed: {sorted(agg.get('formed', []))}")
print(f"   Spectator: {sorted(agg.get('spectator', []))}")

print(f"\n5. DETECTED REACTIONS:")
detection = detect_reaction_types(rxn)
for i, m in enumerate(detection.matches[:5]):
    print(f"   {i+1}. {m.reaction_type}")
    print(f"      Slots: {m.matched_slots}/{m.required_slots}")
    print(f"      Evidence: {m.slot_evidence}")

print("\n" + "=" * 70)
print("ANALYSIS:")
print("  Expected: C_N_Coupling (HeteroAr-Br + amine → HeteroAr-NHR)")
print("  Problem: Ar-Ar is SPECTATOR (biphenyl already present)")
print("           But Arylation_Ar_H matches using it as PRODUCT")
print("=" * 70)
