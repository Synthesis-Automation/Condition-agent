from typing import Dict, Any
import re
def _guess_lg(s: str) -> str:
    s=s.lower()
    if "os(=o)(=o)c(f)(f)f" in s or "otf" in s: return "OTf"
    if "i" in s: return "I"
    if "br" in s: return "Br"
    if "cl" in s: return "Cl"
    return "UNK"
def _nuc_class(s: str) -> str:
    t=s.lower()
    if "indole" in t or "[nh]" in t and "c1" in t: return "indole"
    if re.search(r"c[^)]*n", t): return "aniline"
    if "n(" in t: return "amine_secondary"
    if "n" in t: return "amine_primary"
    if "o" in t and "n" not in t: return "phenol"
    return "amine"
def featurize(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    LG=_guess_lg(electrophile); nuc=_nuc_class(nucleophile)
    ortho_count=0; para_EWG=any(x in electrophile.lower() for x in ["[n+](=o)[o-]","c#n","c(f)(f)f"])
    return {"LG":LG,"elec_class":"aryl","ortho_count":ortho_count,"para_EWG":para_EWG,"heteroaryl":False,
            "nuc_class":nuc,"n_basicity":"aromatic_primary" if nuc in ("aniline",) else "unknown",
            "steric_alpha":"low","bin":f"LG:{LG}|NUC:{nuc}"}
