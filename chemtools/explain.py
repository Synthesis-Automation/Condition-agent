from typing import Dict, Any
def for_pack(pack: Dict[str, Any], features: Dict[str, Any]) -> Dict[str, Any]:
    bin_str = features.get("bin") or f"LG:{features.get('LG','?')}|NUC:{features.get('nuc_class','?')}"
    reasons = [f"Matches substrate bin {bin_str}.", "Core/tails frequently co-used in Ullmann literature for this bin."]
    precedents = pack.get("precedents", [])
    return {"reasons": reasons, "precedents": precedents[:3]}
