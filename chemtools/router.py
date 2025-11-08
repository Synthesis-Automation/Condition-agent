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
    RCM_FAMILIES_CANONICAL,
    ULLMANN_SPECIFIC_CANONICAL,
    apply_catalyst_override,
    canonical_family_label,
    resolve_reaction_family as _resolve_reaction_family,
)
from .util import functional_groups as _functional_groups
from .util.rdkit_helpers import parse_smiles
from .smiles import normalize_reaction as _normalize_reaction


_ROUTER_GROUPS = (
    "acid",
    "acyl_halide",
    "alcohol",
    "aldehyde",
    "alkene",
    "alkoxide",
    "alkyl_halide",
    "alpha_beta_unsaturated",
    "aryl_halide",
    "borane",
    "boron",
    "carbonyl",
    "conjugated_diene",
    "cyanide",
    "diene",
    "ester",
    "grignard",
    "iodide",
    "ketone",
    "nucleophile_n",
    "nucleophile_o",
    "nucleophile_s",
    "organolithium",
    "organozinc",
    "phenol",
    "terminal_alkene",
    "terminal_alkyne",
    "triflate",
    "vinyl_halide",
)


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
    token_hits = _functional_groups.detect_any(reactants, group_subset=_ROUTER_GROUPS)
    result: Dict[str, bool] = {}
    for legacy in _ROUTER_GROUPS:
        token = _functional_groups._resolve_group_token(legacy)  # type: ignore[attr-defined]
        if token is None:
            result[legacy] = False
        else:
            result[legacy] = bool(token_hits.get(token))
    return result


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


def _detect_radical_initiator(reactants: List[str]) -> Optional[str]:
    """Detect radical initiators and halogen sources for radical reactions."""
    reactant_text = " ".join(reactants).lower()
    
    # Radical initiators
    if "aibn" in reactant_text or "azobisisobutyronitrile" in reactant_text:
        return "AIBN"
    if "bpo" in reactant_text or "benzoyl peroxide" in reactant_text:
        return "BPO"
    
    # Halogen sources for radical halogenation
    if "ccl4" in reactant_text or "carbon tetrachloride" in reactant_text:
        return "CCl4"
    if "cbr4" in reactant_text or "carbon tetrabromide" in reactant_text:
        return "CBr4"
    if "nbs" in reactant_text or "n-bromosuccinimide" in reactant_text:
        return "NBS"
    
    # Check SMILES patterns for halogen sources
    for r in reactants:
        r_lower = r.lower()
        # CCl4, CBr4 patterns
        if "c(cl)(cl)(cl)cl" in r_lower or "ccl4" in r_lower:
            return "CCl4"
        if "c(br)(br)(br)br" in r_lower or "cbr4" in r_lower:
            return "CBr4"
        # Br2, Cl2
        if r in ("BrBr", "[Br][Br]", "ClCl", "[Cl][Cl]"):
            return "Br2/Cl2"
    
    return None

