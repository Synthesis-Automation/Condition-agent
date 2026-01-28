"""Debug the full featurization pipeline to see what's happening."""

import sys
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.featurizers.unified import featurize_reaction

# Your reaction
reaction_smiles = "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1"

# Featurize with extended options
options = {
    "confirm_coupling_products": True,
    "detailed": True,
    "motif_site_filter": "substituent",
}

payload = featurize_reaction(reaction_smiles, options=options)

print("=" * 80)
print("FULL PAYLOAD (formatted)")
print("=" * 80)
print(json.dumps(payload, indent=2, sort_keys=False))
