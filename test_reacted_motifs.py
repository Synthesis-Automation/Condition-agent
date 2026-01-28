"""Check what's in reacted_motifs"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction

result = featurize_reaction(
    'c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1',
    options={'detailed': True, 'confirm_coupling_products': True}
)

print("Top-level keys:", list(result.keys()))

extended = result.get('extended', {})
aggregates_ext = extended.get('aggregates', {})

print("\nExtended aggregates:")
print(f"  Reacted: {aggregates_ext.get('reacted_motifs')}")

# Also check top level
print("\nTop-level aggregates:")
core_agg = result.get('aggregates', {})
if core_agg:
    print(f"  Reacted: {core_agg.get('reacted_motifs')}")
else:
    print("  No top-level aggregates")


