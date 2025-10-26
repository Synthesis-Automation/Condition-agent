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
        "carbonyl": Chem.MolFromSmarts("[CX3](=O)[C,H]"),  # Ketone or aldehyde
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


def detect_family(reactants: List[str]) -> Dict[str, Any]:
    h = _rule_hits(reactants)
    fam: Optional[str] = None
    conf = 0.3
    
    # Helper for checking reactant text
    rs = " ".join(reactants).lower()

    # Determine family based on prioritized deterministic rules
    is_aryl_or_vinyl_electrophile = h.get("aryl_halide") or h.get("vinyl_halide") or h.get("triflate")
    
    # PRIORITY 1: Grignard/organometallic addition to carbonyl (HIGHEST PRIORITY)
    # This must come BEFORE any C-C coupling to prioritize carbonyl addition
    if h.get("carbonyl") and (h.get("grignard") or h.get("organolithium") or h.get("organozinc")):
        if h.get("grignard"):
            fam, conf = "grignard_addition", 0.90
        elif h.get("organozinc"):
            fam, conf = "organozinc_addition", 0.90
        else:
            fam, conf = "organolithium_addition", 0.90
    
    # PRIORITY 2: C-C Couplings (if not carbonyl addition)
    # C-C Suzuki Coupling (highest confidence C-C coupling)
    elif h.get("aryl_halide") and h.get("boron"):
        fam, conf = "suzuki_miyaura", 0.9
    
    # C-C Sonogashira Coupling
    elif is_aryl_or_vinyl_electrophile and h.get("terminal_alkyne"):
        fam, conf = "sonogashira", 0.85
    
    # Kumada Coupling (aryl halide + Grignard, no carbonyl already checked above)
    elif is_aryl_or_vinyl_electrophile and h.get("grignard"):
        fam, conf = "kumada_coupling", 0.85
    
    # Negishi Coupling (aryl halide + organozinc, no carbonyl already checked above)
    elif is_aryl_or_vinyl_electrophile and h.get("organozinc"):
        fam, conf = "negishi_coupling", 0.85
    
    # Heck Coupling (aryl halide + alkene, no boron)
    elif is_aryl_or_vinyl_electrophile and h.get("alkene") and not h.get("boron"):
        fam, conf = "heck_coupling", 0.80
    
    # PRIORITY 3: C-N/C-O/C-S Couplings (heteroatom couplings)
    # C-N Coupling (unified - metal preference handled via constraints)
    elif is_aryl_or_vinyl_electrophile and h.get("nucleophile_n"):
        fam, conf = "cn_coupling", 0.9 if h.get("aryl_halide") else 0.8
    
    # C-O Coupling (Ullmann-type etherification)
    elif is_aryl_or_vinyl_electrophile and h.get("nucleophile_o"):
        fam, conf = "co_coupling", 0.85 if h.get("aryl_halide") else 0.75
    
    # C-S Coupling (Ullmann-type thioetherification)
    elif is_aryl_or_vinyl_electrophile and h.get("nucleophile_s"):
        fam, conf = "cs_coupling", 0.85 if h.get("aryl_halide") else 0.75
    
    # PRIORITY 4: Amide/Ester formation
    # Amide formation (higher priority than esterification)
    elif h.get("acid") and h.get("nucleophile_n"):
        fam, conf = "amide_coupling", 0.8
    
    # Esterification (acid + alcohol, but not amide)
    elif h.get("acid") and h.get("alcohol"):
        fam, conf = "esterification", 0.85
    
    # PRIORITY 5: SN2 reactions and substitutions
    # Hydroboration (alkene + borane)
    elif h.get("alkene") and h.get("borane"):
        fam, conf = "hydroboration", 0.85
    
    # SN2 with cyanide (highest priority SN2)
    elif h.get("alkyl_halide") and h.get("cyanide"):
        fam, conf = "nitrile_formation", 0.90
    
    # Finkelstein (halide exchange)
    elif h.get("alkyl_halide") and h.get("iodide"):
        fam, conf = "finkelstein", 0.85
    
    # Williamson ether synthesis
    elif h.get("alkyl_halide") and (h.get("alkoxide") or h.get("alcohol")):
        fam, conf = "williamson_ether", 0.85
    
    # PRIORITY 6: Condensation reactions (LOWER PRIORITY)
    elif h.get("ester") and rs.count("c(=o)o") >= 2:  # Two esters (Claisen)
        fam, conf = "claisen_condensation", 0.70
    
    elif h.get("alpha_beta_unsaturated"):  # α,β-unsaturated carbonyl - Michael
        fam, conf = "michael_addition", 0.65
    
    # PRIORITY 7: Cycloaddition
    elif h.get("conjugated_diene") and h.get("alkene"):
        fam, conf = "diels_alder", 0.85

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
        or fam_rule in CN_FAMILIES_CANONICAL
    )
    fam_rule = apply_catalyst_override(fam_rule, agent_metals, is_cn_coupling=has_cn_signature)
    if fam_rule in CN_FAMILIES_CANONICAL:
        conf_rule = max(conf_rule, 0.9 if hits.get("aryl_halide") else 0.85)
    
    # Phase 2 additions - Reagent-based detection (hydrogenation, reduction, oxidation, elimination)
    # Check reactants for common reagents (they appear in reactants, not separate agents)
    # ONLY override if no family detected yet OR low confidence
    
    # Detect hydrogenation (H2 + alkene/alkyne/aryl/carbonyl)
    reducing_agent = _detect_reducing_agent(reactants)
    if reducing_agent == "H2" and (hits.get("alkene") or hits.get("terminal_alkyne") or hits.get("carbonyl")):
        if fam_rule == "Unknown" or conf_rule < 0.85:
            fam_rule, conf_rule = "hydrogenation", 0.90
    
    # Detect carbonyl reduction (NaBH4/LiAlH4 + carbonyl)
    elif reducing_agent in ("NaBH4", "LiAlH4", "BH3", "DIBAL") and hits.get("carbonyl"):
        if fam_rule == "Unknown" or conf_rule < 0.85:
            fam_rule, conf_rule = "carbonyl_reduction", 0.88
    
    # Detect alcohol oxidation (oxidant + alcohol)
    oxidizing_agent = _detect_oxidizing_agent(reactants)
    if oxidizing_agent in ("PCC", "KMnO4", "CrO3", "Swern", "MnO2", "Dess-Martin") and hits.get("alcohol"):
        if fam_rule == "Unknown" or conf_rule < 0.80:
            fam_rule, conf_rule = "alcohol_oxidation", 0.85
    
    # Detect epoxidation (mCPBA/H2O2 + alkene)
    elif oxidizing_agent in ("mCPBA", "H2O2") and hits.get("alkene"):
        if fam_rule == "Unknown" or conf_rule < 0.80:
            fam_rule, conf_rule = "epoxidation", 0.85
    
    # Detect E2 elimination (strong base + alkyl halide)
    if _detect_strong_base(reactants) and hits.get("alkyl_halide") and not hits.get("cyanide"):
        if fam_rule == "Unknown" or conf_rule < 0.75:
            fam_rule, conf_rule = "e2_elimination", 0.80

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
        fam_rxn_canonical = canonical_family_label(fam_rxn)
        if fam_rxn_canonical:
            fam_rxn = fam_rxn_canonical
        fam_rxn = apply_catalyst_override(
            fam_rxn,
            agent_metals,
            is_cn_coupling=has_cn_signature or fam_rxn in CN_FAMILIES_CANONICAL,
        )
        if fam_rxn in CN_FAMILIES_CANONICAL and (conf_rxn or 0) < 0.85:
            conf_rxn = 0.85

    # Choose final family with intelligent preference logic
    # Prefer rule-based detection for C-O and C-S coupling (more specific)
    # Prefer rxn_insight for other reactions when available and has confidence
    if (fam_rule in CO_FAMILIES_CANONICAL or fam_rule in CS_FAMILIES_CANONICAL) and conf_rule >= 0.75:
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
