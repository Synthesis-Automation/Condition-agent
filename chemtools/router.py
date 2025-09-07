from typing import List, Dict, Any
def _rule_hits(reactants: List[str]) -> Dict[str, bool]:
    rs = ' '.join(reactants).lower()
    return {
        "aryl_halide": any(x in rs for x in ["cl","br","i"]) and "c1" in rs,
        "triflate": "os(=o)(=o)c(f)(f)f" in rs or "otf" in rs,
        "boron": "b(" in rs or "b[" in rs,
        "nucleophile_n": "n" in rs,
        "acid": "c(=o)oh" in rs or "c(=o)o" in rs,
    }
def detect_family(reactants: List[str]) -> Dict[str, Any]:
    h = _rule_hits(reactants); fam="Unknown"; conf=0.3
    if (h["aryl_halide"] or h["triflate"]) and h["nucleophile_n"]: fam, conf = "Ullmann_CN", 0.8
    if h["aryl_halide"] and h["boron"]: fam, conf = "Suzuki_CC", 0.85
    if h["acid"] and h["nucleophile_n"]: fam, conf = "Amide_Coupling", max(conf, 0.7)
    return {"family": fam, "confidence": conf, "hits": h}
