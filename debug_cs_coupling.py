#!/usr/bin/env python
"""Debug: Why doesn't dithiocarbamate + aryl iodide detect as C_S_Coupling?"""

from chemtools.featurizers.formatters.reaction import featurize_reaction
from chemtools.reaction_type_detection import detect_reaction_types

rxn = 'CN(C)C(=S)S.[Na].Clc1ccc(I)cc1>>CN(C)C(=S)Sc1ccc(Cl)cc1'

print("=" * 70)
print("Reaction: Dithiocarbamate + aryl iodide → aryl sulfide")
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
print(f"   Reacted: {sorted(agg.get('reacted_motifs', []))}")
print(f"   Formed: {sorted(agg.get('formed_motifs', []))}")
print(f"   Spectator: {sorted(agg.get('spectator_motifs', []))}")

print(f"\n5. DETECTED REACTIONS:")
detection = detect_reaction_types(rxn)
if detection.matches:
    for i, m in enumerate(detection.matches[:5]):
        print(f"   {i+1}. {m.reaction_type}: {m.slot_evidence}")
else:
    print("   ❌ NO MATCHES")

print("\n" + "=" * 70)
print("EXPECTED: C_S_Coupling (Ar-I + thiol/thiolate → Ar-SR)")
print("=" * 70)
