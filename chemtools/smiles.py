from typing import Dict, Any
import re
from .util.rdkit_helpers import canonical_smiles

def _largest_fragment(smiles: str) -> str:
    parts = smiles.split('.'); return max(parts, key=len)

def normalize(smi: str) -> Dict[str, Any]:
    smi = (smi or '').strip(); fragments = smi.split('.'); largest = _largest_fragment(smi)
    can = canonical_smiles(largest); out = can or largest
    out = out.replace('[Na+]', '').replace('[K+]', '').replace('[Cl-]', 'Cl')
    out = re.sub(r'\[H\]|\[h\]', '', out)
    return {"input": smi, "fragments": fragments, "largest_smiles": largest, "smiles_norm": out}
