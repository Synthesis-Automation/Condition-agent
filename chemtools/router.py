from typing import List, Dict, Any

from .util.rdkit_helpers import rdkit_available, parse_smiles


def _compile_smarts():
    if not rdkit_available():
        return None
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return None
    smarts = {
        # Aryl halide: aromatic carbon bound to Cl/Br/I
        "aryl_halide": Chem.MolFromSmarts("[$(c[Cl,Br,I]),$(c-[Cl,Br,I])]") ,
        # Vinyl halide/triflate (simple patterns)
        "vinyl_halide": Chem.MolFromSmarts("C=C[Cl,Br,I]"),
        "triflate": Chem.MolFromSmarts("OS(=O)(=O)C(F)(F)F"),
        # Boron partners (boronic acids/esters)
        "boron": Chem.MolFromSmarts("[BX3;$(B(O)O),$(B(O)O),$(B(O)O)]"),
        # Terminal alkyne
        "terminal_alkyne": Chem.MolFromSmarts("C#C[H]") or Chem.MolFromSmarts("[C;H]#C"),
        # Carboxylic acid
        "acid": Chem.MolFromSmarts("C(=O)[OH]"),
        # N-nucleophile (amine/anilines, simple)
        "nucleophile_n": Chem.MolFromSmarts("[NX3;H1,H2]"),
        # Phenoxide/alcohol O-H (for etherification heuristic if needed)
        "nucleophile_o": Chem.MolFromSmarts("[OX2H]")
    }
    return smarts


_SMARTS = _compile_smarts()


def _rule_hits(reactants: List[str]) -> Dict[str, bool]:
    # RDKit SMARTS matching when available
    if _SMARTS is not None and rdkit_available():
        try:
            from rdkit import Chem  # type: ignore
        except Exception:
            pass
        mols = [parse_smiles(s) for s in reactants]
        mols = [m for m in mols if m is not None]
        def any_match(key: str) -> bool:
            patt = _SMARTS.get(key)
            if patt is None:
                return False
            for m in mols:
                try:
                    if m.HasSubstructMatch(patt):
                        return True
                except Exception:
                    continue
            return False
        return {
            "aryl_halide": any_match("aryl_halide"),
            "vinyl_halide": any_match("vinyl_halide"),
            "triflate": any_match("triflate"),
            "boron": any_match("boron"),
            "nucleophile_n": any_match("nucleophile_n"),
            "nucleophile_o": any_match("nucleophile_o"),
            "terminal_alkyne": any_match("terminal_alkyne"),
            "acid": any_match("acid"),
        }

    # Text fallback heuristic when RDKit not available
    rs = " ".join(reactants).lower()
    def has(pattern: str) -> bool:
        return pattern in rs
    return {
        "aryl_halide": (has("cl") or has("br") or has(" i")) and ("c1" in rs or "c2" in rs or "c[" in rs),
        "vinyl_halide": (has("c=ccl") or has("c=cbr") or has("c=ci")),
        "triflate": has("os(=o)(=o)c(f)(f)f") or has("otf"),
        "boron": has("b(") or has("b[") or has("b(o)o"),
        "nucleophile_n": has("n") or has("nh"),
        "nucleophile_o": has("o") or has("oh"),
        "terminal_alkyne": has("c#c") or has("c#cc"),
        "acid": has("c(=o)oh") or has("c(=o)o") or has("oc(=o)"),
    }


def detect_family(reactants: List[str]) -> Dict[str, Any]:
    h = _rule_hits(reactants)
    fam = "Unknown"
    conf = 0.3

    # Determine family based on prioritized deterministic rules
    is_aryl_or_vinyl_electrophile = h.get("aryl_halide") or h.get("vinyl_halide") or h.get("triflate")
    if is_aryl_or_vinyl_electrophile and h.get("nucleophile_n"):
        fam, conf = "Ullmann_CN", 0.9 if h.get("aryl_halide") else 0.8
    if h.get("aryl_halide") and h.get("boron"):
        fam, conf = "Suzuki_CC", max(conf, 0.9)
    if is_aryl_or_vinyl_electrophile and h.get("terminal_alkyne"):
        fam, conf = "Sonogashira_CC", max(conf, 0.85)
    if h.get("acid") and h.get("nucleophile_n"):
        fam, conf = "Amide_Coupling", max(conf, 0.8)
    # Optional: etherification (Ullmann O) – not primary target but can hint
    if is_aryl_or_vinyl_electrophile and h.get("nucleophile_o") and fam == "Unknown":
        fam, conf = "Ullmann_O", 0.75

    return {"family": fam, "confidence": float(conf), "hits": h}
