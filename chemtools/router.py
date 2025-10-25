from typing import List, Dict, Any, Optional, Set

import re

from .taxonomy import load_registry
from .taxonomy.registry import TaxonomyRegistry
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
        "nucleophile_s": Chem.MolFromSmarts("[SX2H]")
    }
    return smarts


_SMARTS = _compile_smarts()

_TAXONOMY_CACHE: Optional[TaxonomyRegistry] = None
_TAXONOMY_LOAD_FAILED = False

_FAMILY_ALIAS_OVERRIDES: Dict[str, str] = {
    "C_N_Coupling": "cn_coupling",
    "Buchwald_CN": "buchwald_hartwig_c_n",
    "Ullmann_CN": "ullmann_cn",
    "Chan_Lam_CN": "chan_lam_cn",
    "C_O_Coupling": "co_coupling",
    "Ullmann_CO": "ullmann_ether",
    "C_S_Coupling": "cs_coupling",
    "Suzuki_CC": "suzuki_miyaura",
    "Negishi": "negishi",
    "Sonogashira_CC": "sonogashira",
    "Stille": "stille",
    "Heck": "heck",
    "Amide_Coupling": "amide_coupling",
}

_CN_FAMILIES_CANONICAL = {"cn_coupling", "buchwald_hartwig_c_n", "ullmann_cn", "chan_lam_cn"}
_CO_FAMILIES_CANONICAL = {"co_coupling", "ullmann_ether"}
_CS_FAMILIES_CANONICAL = {"cs_coupling"}
_AMIDE_FAMILIES_CANONICAL = {"amide_coupling"}
_SONOGASHIRA_FAMILIES_CANONICAL = {"sonogashira"}
_SUZUKI_FAMILIES_CANONICAL = {"suzuki_miyaura"}
_ULLMANN_SPECIFIC_CANONICAL = {"ullmann_cn"}
_BUCHWALD_SPECIFIC_CANONICAL = {"buchwald_hartwig_c_n"}


def _get_taxonomy_registry() -> Optional[TaxonomyRegistry]:
    global _TAXONOMY_CACHE, _TAXONOMY_LOAD_FAILED
    if _TAXONOMY_CACHE is not None:
        return _TAXONOMY_CACHE
    if _TAXONOMY_LOAD_FAILED:
        return None
    try:
        _TAXONOMY_CACHE = load_registry()
        return _TAXONOMY_CACHE
    except Exception:
        _TAXONOMY_LOAD_FAILED = True
        return None


def _slugify_family(value: str) -> str:
    return re.sub(r"[^0-9a-z]+", "_", value.lower()).strip("_")


def _canonical_family_label(family: Optional[str]) -> Optional[str]:
    if not family:
        return None
    family = family.strip()
    if not family:
        return None

    registry = _get_taxonomy_registry()
    if registry:
        if family in registry.reaction_types:
            return family
        alias = registry.resolve_alias(family)
        if alias and alias.entity_type == "reaction_type":
            return alias.entity_id
        alias = registry.resolve_alias(family.lower())
        if alias and alias.entity_type == "reaction_type":
            return alias.entity_id
        slug = _slugify_family(family)
        alias = registry.resolve_alias(slug)
        if alias and alias.entity_type == "reaction_type":
            return alias.entity_id
    return _FAMILY_ALIAS_OVERRIDES.get(family)


def resolve_reaction_family(family: Optional[str]) -> Optional[str]:
    """Resolve arbitrary family label/alias to canonical taxonomy ID (or None)."""
    return _canonical_family_label(family)

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


def _apply_catalyst_override(
    family: str,
    metals: Set[str],
    *,
    is_cn_coupling: bool,
) -> str:
    """
    Apply catalyst-based family override for C-N coupling reactions using canonical IDs.
    """
    canonical = family

    if not metals or not is_cn_coupling:
        return canonical

    if "Pd" in metals and canonical not in _BUCHWALD_SPECIFIC_CANONICAL:
        override = _canonical_family_label("Buchwald_CN")
        if override:
            return override
    if "Cu" in metals and canonical not in _ULLMANN_SPECIFIC_CANONICAL:
        override = _canonical_family_label("Ullmann_CN")
        if override:
            return override

    return canonical


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
        }
        # Combine conservatively (logical OR)
        return {k: bool(rdkit_hits.get(k) or text_hits.get(k)) for k in text_hits.keys()}

    return text_hits


def detect_family(reactants: List[str]) -> Dict[str, Any]:
    h = _rule_hits(reactants)
    fam: Optional[str] = None
    conf = 0.3

    # Determine family based on prioritized deterministic rules
    is_aryl_or_vinyl_electrophile = h.get("aryl_halide") or h.get("vinyl_halide") or h.get("triflate")
    
    # C-N Coupling (unified - metal preference handled via constraints)
    if is_aryl_or_vinyl_electrophile and h.get("nucleophile_n"):
        fam, conf = "cn_coupling", 0.9 if h.get("aryl_halide") else 0.8
    
    # C-O Coupling (Ullmann-type etherification)
    if is_aryl_or_vinyl_electrophile and h.get("nucleophile_o"):
        fam, conf = "co_coupling", 0.85 if h.get("aryl_halide") else 0.75
    
    # C-S Coupling (Ullmann-type thioetherification)
    if is_aryl_or_vinyl_electrophile and h.get("nucleophile_s"):
        fam, conf = "cs_coupling", 0.85 if h.get("aryl_halide") else 0.75
    
    # C-C Suzuki Coupling (higher priority than C-O/C-S)
    if h.get("aryl_halide") and h.get("boron"):
        fam, conf = "suzuki_miyaura", max(conf, 0.9)
    
    # C-C Sonogashira Coupling
    if is_aryl_or_vinyl_electrophile and h.get("terminal_alkyne"):
        fam, conf = "sonogashira", max(conf, 0.85)
    
    # Amide formation
    if h.get("acid") and h.get("nucleophile_n"):
        fam, conf = "amide_coupling", max(conf, 0.8)

    canonical_family = fam or "Unknown"

    return {
        "family": canonical_family,
        "confidence": float(conf),
        "hits": h,
    }


def detect_family_from_reaction(reaction_smiles: str, *, use_rxn_insight: bool = True) -> Dict[str, Any]:
    """Detect reaction family from a reaction SMILES string.

    - If ``use_rxn_insight`` and the optional rxn_insight integration is available,
      try it first and map to our family labels; otherwise fall back to rule-based
      detection over reactant SMARTS.
    - Always returns a dict with keys: family, confidence, hits, auto (optional).
    """
    rsmi = (reaction_smiles or "").strip()
    norm = _normalize_reaction(rsmi)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]
    reactants = [s for s in reactants if s]
    agent_metals = _detect_agent_metals(norm.get("agents") or [])

    # Base rule-based detection (deterministic, no network)
    base = detect_family(reactants)
    fam_rule = base.get("family") or "Unknown"
    conf_rule: float = float(base.get("confidence") or 0.3)
    hits = base.get("hits") or _rule_hits(reactants)
    if agent_metals:
        hits = dict(hits)
        hits["catalyst_pd"] = "Pd" in agent_metals
        hits["catalyst_cu"] = "Cu" in agent_metals
        hits["catalyst_ni"] = "Ni" in agent_metals
        hits["catalyst_co"] = "Co" in agent_metals

    has_cn_signature = bool(
        (hits.get("nucleophile_n") and (hits.get("aryl_halide") or hits.get("vinyl_halide") or hits.get("triflate")))
        or fam_rule in _CN_FAMILIES_CANONICAL
    )
    fam_rule = _apply_catalyst_override(fam_rule, agent_metals, is_cn_coupling=has_cn_signature)
    if fam_rule in _CN_FAMILIES_CANONICAL:
        conf_rule = max(conf_rule, 0.9 if hits.get("aryl_halide") else 0.85)

    auto: Optional[Dict[str, Any]] = None
    fam_rxn: Optional[str] = None
    conf_rxn: Optional[float] = None
    if use_rxn_insight:
        try:
            from .reaction_type_detector import detect_reaction_type as _rxn_detect, is_available as _rxn_avail  # type: ignore
        except Exception:
            _rxn_detect = None  # type: ignore
            _rxn_avail = lambda: False  # type: ignore
        if _rxn_detect is not None and _rxn_avail():
            try:
                auto = _rxn_detect(norm.get("normalized") or rsmi)  # type: ignore[misc]
            except Exception:
                auto = None
            if isinstance(auto, dict) and auto.get("success") and auto.get("mapped_family"):
                fam_rxn = str(auto.get("mapped_family") or "Unknown")
                aconf = auto.get("confidence")
                try:
                    conf_rxn = float(aconf) if aconf is not None else None
                except Exception:
                    conf_rxn = None

    if fam_rxn:
        fam_rxn_canonical = _canonical_family_label(fam_rxn)
        if fam_rxn_canonical:
            fam_rxn = fam_rxn_canonical
        fam_rxn = _apply_catalyst_override(
            fam_rxn,
            agent_metals,
            is_cn_coupling=has_cn_signature or fam_rxn in _CN_FAMILIES_CANONICAL,
        )
        if fam_rxn in _CN_FAMILIES_CANONICAL and (conf_rxn or 0) < 0.85:
            conf_rxn = 0.85

    # Choose final family with intelligent preference logic
    # Prefer rule-based detection for C-O and C-S coupling (more specific)
    # Prefer rxn_insight for other reactions when available and has confidence
    if (fam_rule in _CO_FAMILIES_CANONICAL or fam_rule in _CS_FAMILIES_CANONICAL) and conf_rule >= 0.75:
        # Rule-based C-O/C-S detection is reliable, use it
        fam_final = fam_rule
        conf_final = conf_rule
    elif fam_rxn and conf_rxn is not None:
        # RXN insight has a result with confidence
        fam_final = fam_rxn
        conf_final = conf_rxn
    else:
        # Fall back to rule-based
        fam_final = fam_rule
        conf_final = conf_rule
    
    agreement = (fam_rxn is not None) and (str(fam_rxn) == str(fam_rule))

    out: Dict[str, Any] = {
        "family": fam_final,
        "confidence": float(conf_final),
        "hits": hits,
        "agreement": bool(agreement),
        "status": ("consistent" if agreement else ("conflict" if fam_rxn else "rule_only")),
    }
    if auto is not None:
        out["auto"] = auto
    # Attach both results for comparison
    out["rule"] = {
        "family": fam_rule,
        "confidence": conf_rule,
        "hits": hits,
    }
    if fam_rxn or auto is not None:
        out["rxn"] = {
            "family": fam_rxn,
            "confidence": conf_rxn,
            **({} if auto is None else {
                "available": bool(auto.get("available")),
                "success": bool(auto.get("success")),
                "rxn_class": auto.get("rxn_class"),
                "rxn_name": auto.get("rxn_name"),
            }),
        }
    if agent_metals:
        out["catalysts"] = {
            "metals": sorted(agent_metals),
            "source": "agents",
        }
    return out
