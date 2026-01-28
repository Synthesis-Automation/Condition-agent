"""
Investigate why Suzuki reaction is classified as Arylation_Ar_H
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction
import json

result = featurize_reaction(
    'c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1',
    options={'detailed': True, 'confirm_coupling_products': True}
)

print("=" * 70)
print("REACTION TYPE DETECTION ANALYSIS")
print("=" * 70)

print("\n=== MAIN REACTION TYPE (from core) ===")
print(f"reaction_type: {result.get('reaction_type')}")
print(f"confidence: {result.get('confidence')}")
print(f"reaction_key: {result.get('reaction_key')}")

extended = result.get('extended', {})

print("\n=== DETECTION MATCHES (from extended) ===")
detection = extended.get('detection', {})
matches = detection.get('matches', [])
print(f"Total matches: {len(matches)}")
if matches:
    for i, match in enumerate(matches[:5], 1):
        print(f"  {i}. {match.get('name')} - confidence: {match.get('confidence'):.4f}")
        slot_evidence = match.get('slot_evidence', {})
        if slot_evidence:
            print(f"     Slots: {list(slot_evidence.keys())}")

print("\n=== ROLE CLASSIFICATION (from extended) ===")
role_class = extended.get('role_classification', {})
roles = role_class.get('reactants', {})
print(f"reaction_type: {roles.get('reaction_type')}")
print(f"confidence: {roles.get('confidence')}")
print(f"detection_method: {roles.get('detection_method')}")
print(f"num_reactants: {roles.get('num_reactants')}")

reactants_list = roles.get('reactants', [])
print(f"\nReactant roles:")
for r in reactants_list:
    print(f"  - {r.get('category')} ({r.get('role')})")

print("\n=== AGGREGATES (functional groups) ===")
aggregates = extended.get('aggregates', {})
print(f"reacted_motifs: {aggregates.get('reacted_motifs')}")
print(f"formed_motifs: {aggregates.get('formed_motifs')}")
print(f"motif_ids: {aggregates.get('motif_ids')}")

print("\n" + "=" * 70)
print("ANALYSIS")
print("=" * 70)
print("""
The reaction has:
- Ar-Br (aryl bromide)  
- Ar-B(OH)2 (boronic acid)
- Forms Ar-Ar (biaryl)

This is a classic Suzuki-Miyaura coupling!

Possible reasons for misclassification:
1. Detection system prioritizes different rules
2. "Arylation_Ar_H" may be matching before Suzuki
3. Slot-based detection may be misconfigured
""")
