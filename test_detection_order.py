"""Quick test of detection with coupling confirmation"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[0]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.detection.core import detect_reaction_types

det = detect_reaction_types(
    'c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1',
    confirm_coupling_products=True
)

print('Detection matches with confirm_coupling_products=True:')
for i, m in enumerate(det.matches[:5], 1):
    print(f'  {i}. {m.reaction_type} (confidence: {m.confidence})')
    print(f'     Slots: {list(m.slot_evidence.keys())}')
