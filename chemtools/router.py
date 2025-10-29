from typing import Any, Dict, List, Optional, Set

import re

from .analysis.reactions import (
    AMIDE_FAMILIES_CANONICAL,
    BUCHWALD_SPECIFIC_CANONICAL,
    CN_FAMILIES_CANONICAL,
    CO_FAMILIES_CANONICAL,
    CS_FAMILIES_CANONICAL,
    SONOGASHIRA_FAMILIES_CANONICAL,
    SUZUKI_FAMILIES_CANONICAL,
    ULLMANN_SPECIFIC_CANONICAL,
    apply_catalyst_override,
    canonical_family_label,
    resolve_reaction_family as _resolve_reaction_family,
)
from .util.rdkit_helpers import rdkit_available, parse_smiles
from .smiles import normalize_reaction as _normalize_reaction


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
        # Phenoxide/alcohol O-H (for C-O coupling)
        "nucleophile_o": Chem.MolFromSmarts("[OX2H]"),
        # Thiol S-H (for C-S coupling)
        "nucleophile_s": Chem.MolFromSmarts("[SX2H]"),
        
        # Phase 2 additions - carbonyl compounds
        "carbonyl": Chem.MolFromSmarts("[CX3]=O"),  # Ketone, aldehyde, or ester carbonyl
        "aldehyde": Chem.MolFromSmarts("[CX3H](=O)"),  # Aldehyde specifically
        "ketone": Chem.MolFromSmarts("[CX3](=O)[C]"),  # Ketone specifically
        "ester": Chem.MolFromSmarts("[CX3](=O)[OX2][C,H]"),  # Ester
        
        # Phase 2 additions - nucleophiles and organometallics
        "alcohol": Chem.MolFromSmarts("[CX4][OX2H]"),  # Aliphatic alcohol
        "grignard": Chem.MolFromSmarts("[C,c][Mg][Br,Cl,I]"),  # Grignard reagent
        "organozinc": Chem.MolFromSmarts("[C,c][Zn][Br,Cl,I]"),  # Organozinc (Negishi)
        "organolithium": Chem.MolFromSmarts("[C,c][Li]"),  # Organolithium
        "cyanide": Chem.MolFromSmarts("[C-]#N"),  # Cyanide anion
        "iodide": Chem.MolFromSmarts("[I-]"),  # Iodide anion
        "alkoxide": Chem.MolFromSmarts("[O-][C,H]"),  # Alkoxide anion
        
        # Phase 2 additions - alkyl halides and alkenes
        "alkyl_halide": Chem.MolFromSmarts("[CX4][Cl,Br,I]"),  # Alkyl halide (not aromatic)
        "alkene": Chem.MolFromSmarts("C=C"),  # Simple alkene
        "conjugated_diene": Chem.MolFromSmarts("C=C-C=C"),  # Conjugated diene
        "alpha_beta_unsaturated": Chem.MolFromSmarts("C=C-[CX3]=O"),  # α,β-unsaturated carbonyl
        
        # Phase 2 additions - boron reagents
        "borane": Chem.MolFromSmarts("[BH3,BH2,BH]"),  # Borane (BH3, B2H6, etc.)
    }
    return smarts


_SMARTS = _compile_smarts()


def resolve_reaction_family(family: Optional[str]) -> Optional[str]:
    """Backwards-compatible wrapper around the analysis resolver."""
    return _resolve_reaction_family(family)

_METAL_ATOMIC_NUMBERS = {
    29: "Cu",
    46: "Pd",
    28: "Ni",
    27: "Co",
}

_METAL_NAME_TOKENS = {
    "palladium": "Pd",
    "copper": "Cu",
    "nickel": "Ni",
    "cobalt": "Co",
}

_METAL_REGEXES = {
    "Pd": re.compile(r"(?<![A-Za-z])[Pp][dD](?![a-z])"),
    "Cu": re.compile(r"(?<![A-Za-z])[Cc][uU](?![a-z])"),
    "Ni": re.compile(r"(?<![A-Za-z])[Nn][iI](?![a-z])"),
    "Co": re.compile(r"(?<![A-Za-z])[Cc][oO](?![a-z])"),
}


def _detect_agent_metals(agents: List[Dict[str, Any]]) -> Set[str]:
    """Detect metal catalysts from normalized agent block."""

    metals: Set[str] = set()
    seen_smiles: Set[str] = set()
    for agent in agents or []:
        snippets: List[str] = []
        for key in ("smiles_norm", "largest_smiles", "input"):
            value = agent.get(key)
            if not value:
                continue
            text = str(value).strip()
            if not text or text in seen_smiles:
                continue
            seen_smiles.add(text)
            snippets.append(text)

        for snippet in snippets:
            mol = parse_smiles(snippet)
            if mol is not None:
                try:
                    for atom in mol.GetAtoms():
                        label = _METAL_ATOMIC_NUMBERS.get(atom.GetAtomicNum())
                        if label:
                            metals.add(label)
                except Exception:
                    pass

        joined = " ".join(snippets)
        if joined:
            lower = joined.lower()
            for token, label in _METAL_NAME_TOKENS.items():
                if token in lower:
                    metals.add(label)
            for label, pattern in _METAL_REGEXES.items():
                if pattern.search(joined):
                    metals.add(label)

    return metals


def _rule_hits(reactants: List[str]) -> Dict[str, bool]:
    # Text fallback heuristic (used when RDKit unavailable, or to augment matches when some tokens fail to parse)
    rs = " ".join(reactants).lower()
    def has(pattern: str) -> bool:
        return pattern in rs
    text_hits = {
        "aryl_halide": (has("cl") or has("br") or has(" i")) and ("c1" in rs or "c2" in rs or "c[" in rs),
        "vinyl_halide": (has("c=ccl") or has("c=cbr") or has("c=ci")),
        "triflate": has("os(=o)(=o)c(f)(f)f") or has("otf"),
        # Accept both common orders: 'b(o)o' and 'ob(o)' (RDKit canonicalization may reorder)
        "boron": has("b(") or has("b[") or has("b(o)o") or has("ob(o)"),
        "nucleophile_n": has("n") or has("nh"),
        "nucleophile_o": has("o") or has("oh"),
        "nucleophile_s": has("s") or has("sh"),
        "terminal_alkyne": has("c#c") or has("c#cc"),
        "acid": has("c(=o)oh") or has("c(=o)o") or has("oc(=o)"),
        
        # Phase 2 additions
        "carbonyl": has("c(=o)") or has("c=o"),
        "aldehyde": has("c(=o)") and (has("c(=o)h") or has("ch=o")),
        "ketone": has("c(=o)c") or has("cc(=o)c"),
        "ester": has("c(=o)o") and not has("c(=o)oh"),
        "alcohol": has("co") or has("oh"),
        "grignard": has("[mg]") or has("mgbr") or has("mgcl"),
        "organozinc": has("[zn]") or has("znbr") or has("zncl") or has("rzn"),
        "organolithium": has("[li]") or has("li"),
        "cyanide": has("[c-]#n") or has("cn-") or has("nacn") or has("kcn"),
        "iodide": has("[i-]") or has("nai") or has("ki"),
        "alkoxide": has("[o-]") or has("naome") or has("kome") or has("naoh") or has("koh"),
        "alkyl_halide": (has("cl") or has("br") or has("i")) and has("c"),
        "alkene": has("c=c"),
        "conjugated_diene": has("c=c") and rs.count("c=c") >= 2,
        "alpha_beta_unsaturated": has("c=c") and has("c(=o)"),
        "borane": has("bh3") or has("b2h6") or has("bh2") or has("[bh"),
    }

    # RDKit SMARTS matching when available; OR with text_hits to be robust to parse failures
    if _SMARTS is not None and rdkit_available():
        try:
            from rdkit import Chem  # type: ignore
        except Exception:
            Chem = None  # type: ignore
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
        # RDKit SMARTS matching when available; OR with text_hits to be robust to parse failures
        rdkit_hits = {
            "aryl_halide": any_match("aryl_halide"),
            "vinyl_halide": any_match("vinyl_halide"),
            "triflate": any_match("triflate"),
            "boron": any_match("boron"),
            "nucleophile_n": any_match("nucleophile_n"),
            "nucleophile_o": any_match("nucleophile_o"),
            "nucleophile_s": any_match("nucleophile_s"),
            "terminal_alkyne": any_match("terminal_alkyne"),
            "acid": any_match("acid"),
            # Phase 2 additions
            "carbonyl": any_match("carbonyl"),
            "aldehyde": any_match("aldehyde"),
            "ketone": any_match("ketone"),
            "ester": any_match("ester"),
            "alcohol": any_match("alcohol"),
            "grignard": any_match("grignard"),
            "organozinc": any_match("organozinc"),
            "organolithium": any_match("organolithium"),
            "cyanide": any_match("cyanide"),
            "iodide": any_match("iodide"),
            "alkoxide": any_match("alkoxide"),
            "alkyl_halide": any_match("alkyl_halide"),
            "alkene": any_match("alkene"),
            "conjugated_diene": any_match("conjugated_diene"),
            "alpha_beta_unsaturated": any_match("alpha_beta_unsaturated"),
            "borane": any_match("borane"),
        }
        # Combine conservatively (logical OR)
        return {k: bool(rdkit_hits.get(k) or text_hits.get(k)) for k in text_hits.keys()}

    return text_hits


def _detect_reducing_agent(reactants: List[str]) -> Optional[str]:
    """Detect common reducing agents in reactants. Returns agent type or None."""
    reactant_text = " ".join(reactants).lower()
    
    # Check for hydrogen gas (molecular hydrogen)
    if "[h][h]" in reactant_text or "hh" in reactant_text:
        return "H2"
    # Check for sodium borohydride
    if "nabh4" in reactant_text or "nabh(oac)3" in reactant_text or "[bh4-]" in reactant_text:
        return "NaBH4"
    # Check for lithium aluminum hydride
    if "lialh4" in reactant_text or "[alh4-]" in reactant_text:
        return "LiAlH4"
    # Check for borane
    if "bh3" in reactant_text or "b2h6" in reactant_text:
        return "BH3"
    # Check for DIBAL
    if "dibal" in reactant_text or "diisobutyl" in reactant_text:
        return "DIBAL"
    return None


def _detect_oxidizing_agent(reactants: List[str]) -> Optional[str]:
    """Detect common oxidizing agents in reactants. Returns agent type or None."""
    reactant_text = " ".join(reactants).lower()
    
    # Check for chromium-based oxidants
    if "pcc" in reactant_text or "pyridinium" in reactant_text and "chrom" in reactant_text:
        return "PCC"
    if "kmno4" in reactant_text or "permanganate" in reactant_text or "[mno4-]" in reactant_text:
        return "KMnO4"
    if "cro3" in reactant_text or "jones" in reactant_text:
        return "CrO3"
    # Check for Swern oxidation components
    if "swern" in reactant_text or ("dmso" in reactant_text and ("oxalyl" in reactant_text or "cocl" in reactant_text)):
        return "Swern"
    # Check for manganese dioxide
    if "mno2" in reactant_text:
        return "MnO2"
    # Check for mCPBA (meta-chloroperoxybenzoic acid)
    if "mcpba" in reactant_text or "m-cpba" in reactant_text or "chloroperoxy" in reactant_text:
        return "mCPBA"
    # Check for hydrogen peroxide
    if "h2o2" in reactant_text or "ooh" in reactant_text:
        return "H2O2"
    # Dess-Martin periodinane
    if "dess" in reactant_text or "martin" in reactant_text:
        return "Dess-Martin"
    return None


def _detect_strong_base(reactants: List[str]) -> bool:
    """Detect strong bases suitable for E2 elimination or deprotonation."""
    reactant_text = " ".join(reactants).lower()
    
    strong_bases = ["kot-bu", "kotbu", "koh", "naoh", "nah", "lda", "lhmds", "khmds", "dbu", 
                    "naome", "naoet", "[oh-]", "[o-]", "hydroxide", "alkoxide"]
    return any(base in reactant_text for base in strong_bases)


# Old detection functions removed - use chemtools.detect_reaction() instead
