"""
Analysis: Why Suzuki reaction is misclassified as Arylation_Ar_H

REACTION: c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1

Expected: Suzuki_miyaura
Actual: Arylation_Ar_H
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_molecule

print("=" * 80)
print("ROOT CAUSE ANALYSIS")
print("=" * 80)

# Analyze the boronic acid molecule
mol = featurize_molecule("c1cccnc1B(O)O", options={"detailed": False})

print("\nBoronic acid molecule: c1cccnc1B(O)O")
print("\nDetected motifs:")
for m in mol.get('motifs', []):
    print(f"  - {m['id']} (rank: {m['rank']})")

print("\n" + "=" * 80)
print("THE PROBLEM")
print("=" * 80)
print("""
The boronic acid c1cccnc1B(O)O contains:
1. Ar-B(OH)2 motif (the nucleophile for Suzuki)
2. Multiple Ar-H motifs (aromatic C-H bonds)

Detection Rule Matching:
------------------------

Suzuki_miyaura rule:
  - electrophile: @sp2_electrophiles (Ar-Br) ✓
  - nucleophile: @organoboron (Ar-B(OH)2) ✓
  - product: Ar-Ar ✓

Arylation_Ar_H rule:
  - electrophile: @sp2_electrophiles (Ar-Br) ✓
  - substrate: Ar-H (found in boronic acid!) ✓
  - product: Ar-Ar ✓

Both rules match! The detection system returns the first/highest priority match,
which is currently Arylation_Ar_H.

WHY THIS HAPPENS:
-----------------
The boronic acid molecule HAS aromatic C-H bonds, so it legitimately contains
the "Ar-H" motif. The detection system doesn't know that these C-H bonds are
spectators rather than reactive sites.

SOLUTIONS:
----------
1. Improve motif priority - Ar-B(OH)2 should "consume" or mask Ar-H motifs
   on the same molecule
2. Add exclusion rules - If Ar-B(OH)2 is present, don't match Ar-H substrate
3. Reorder reaction rules - Put Suzuki before Arylation_Ar_H in priority
4. Add slot constraints - Require Ar-H substrate to be on a separate molecule
   from the organoboron nucleophile

WORKAROUND:
-----------
The role classification system correctly identifies this as Suzuki_miyaura
because it uses more sophisticated logic that understands the difference
between nucleophiles and substrates.

""")

print("=" * 80)
print("RECOMMENDATION")
print("=" * 80)
print("""
For conditions recommendation, use the ROLE CLASSIFICATION result, not the
main detection result:

From test output:
  Main detection: Arylation_Ar_H (WRONG)
  Role classification: Suzuki_miyaura (CORRECT)

The role classification system correctly identifies this as a Suzuki reaction
by recognizing that Ar-Br is the electrophile and Ar-B(OH)2 is the nucleophile.
""")
