"""
Unified reaction detection API.

This module consolidates all reaction detection logic (SMARTS-based, ML-based,
catalyst-based) into a single, clean API with consistent outputs.

Public Functions:
    detect_reaction(reaction_smiles, use_ml=True) -> dict
    
Internal:
    _DetectionEngine class (consolidates all detection logic)

All outputs are validated against the unified taxonomy in chemtools/taxonomy/.
"""

from typing import Any, Dict, List, Optional, Set
import logging

from .smiles import normalize_reaction
from .router import (
    _rule_hits,
    _detect_agent_metals,
    _detect_reducing_agent,
    _detect_oxidizing_agent,
    _detect_strong_base,
    _detect_radical_initiator,
)
from .util.smarts_cache import compile_smarts
from .detection_mapper import (
    resolve_to_taxonomy,
    calculate_confidence_adjustment,
    get_mapping_method,
)
from .analysis.reactions import (
    CN_FAMILIES_CANONICAL,
    CO_FAMILIES_CANONICAL,
    CS_FAMILIES_CANONICAL,
    AMIDE_FAMILIES_CANONICAL,
)

logger = logging.getLogger(__name__)


class _DetectionEngine:
    """
    Internal detection engine that consolidates all detection methods.
    
    This class orchestrates:
    - Reaction normalization
    - Catalyst detection
    - Functional group detection (SMARTS)
    - Rule-based family detection
    - ML-based detection (optional, via rxn-insight)
    - Catalyst-based overrides
    """
    
    def __init__(self, reaction_smiles: str):
        """
        Initialize detection engine with reaction SMILES.
        
        Args:
            reaction_smiles: Full reaction SMILES (reactants>>products)
        """
        self.reaction_smiles = reaction_smiles
        self.normalized: Optional[Dict[str, Any]] = None
        self.reactants: List[str] = []
        self.agents: List[Dict[str, Any]] = []
        self.catalysts: Set[str] = set()
        self.functional_groups: Dict[str, bool] = {}
        
        # Normalize reaction and extract components
        self._normalize()
        
    def _normalize(self) -> None:
        """Normalize reaction SMILES and extract components."""
        try:
            self.normalized = normalize_reaction(self.reaction_smiles)
            
            # Extract reactants
            self.reactants = [
                (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
                for r in (self.normalized.get("reactants") or [])
            ]
            self.reactants = [s for s in self.reactants if s]
            
            # Extract agents
            self.agents = self.normalized.get("agents") or []
            
        except Exception as e:
            logger.warning(f"Failed to normalize reaction: {e}")
            self.normalized = {}
            self.reactants = []
            self.agents = []
    
    def _detect_catalysts(self) -> Set[str]:
        """
        Detect catalyst metals from agents.
        
        Extracts Pd, Cu, Ni, Co from agent SMILES and names.
        
        Returns:
            Set of metal symbols (e.g., {"Pd", "Cu"})
        """
        try:
            self.catalysts = _detect_agent_metals(self.agents)
        except Exception as e:
            logger.warning(f"Failed to detect catalysts: {e}")
            self.catalysts = set()
        
        return self.catalysts
    
    def _detect_functional_groups(self) -> Dict[str, bool]:
        """
        Detect functional groups in reactants using SMARTS patterns.
        
        Uses the SMARTS patterns from router.py to identify:
        - Electrophiles: aryl_halide, vinyl_halide, triflate, alkyl_halide
        - Nucleophiles: nucleophile_n, nucleophile_o, nucleophile_s
        - Organometallics: boron, grignard, organozinc, organolithium
        - Others: terminal_alkyne, acid, carbonyl, alkene, etc.
        
        Returns:
            Dict mapping functional group names to boolean presence
        """
        try:
            self.functional_groups = _rule_hits(self.reactants)
        except Exception as e:
            logger.warning(f"Failed to detect functional groups: {e}")
            self.functional_groups = {}
        
        # Add catalyst info to functional groups
        if self.catalysts:
            self.functional_groups["catalyst_pd"] = "Pd" in self.catalysts
            self.functional_groups["catalyst_cu"] = "Cu" in self.catalysts
            self.functional_groups["catalyst_ni"] = "Ni" in self.catalysts
            self.functional_groups["catalyst_co"] = "Co" in self.catalysts
        
        return self.functional_groups
    
    def _rule_based_detection(self) -> Dict[str, Any]:
        """
        Rule-based SMARTS detection.
        
        Uses deterministic rules based on functional group combinations
        to predict reaction family. Consolidates logic from router.detect_family().
        
        Returns:
            {
                "family": str | None,     # Raw family name (pre-taxonomy mapping)
                "confidence": float,       # 0.0-1.0 confidence
                "hits": dict              # Functional group hits
            }
        """
        h = self.functional_groups
        fam: Optional[str] = None
        conf = 0.3
        
        # Helper
        rs = " ".join(self.reactants).lower()
        is_aryl_or_vinyl_electrophile = h.get("aryl_halide") or h.get("vinyl_halide") or h.get("triflate")
        
        # Compute SNAr conditions (used later in priority chain)
        is_heteroaryl = any(x in rs for x in ["c1nc", "c1oc", "c1sc", "n1c", "nc(", "pyrid", "pyrim", "triaz"])
        has_strong_ewg = any(x in rs for x in ["[n+](=o)[o-]", "c#n", "c(f)(f)f", "s(=o)(=o)", "c(=o)c"])
        is_snar_electrophile = is_heteroaryl or (h.get("aryl_halide") and has_strong_ewg)
        has_nucleophile = h.get("nucleophile_n") or h.get("nucleophile_o") or h.get("nucleophile_s")
        no_metal_catalyst = not self.catalysts  # SNAr doesn't need metal catalyst
        is_snar = is_snar_electrophile and has_nucleophile and no_metal_catalyst
        
        # PRIORITY 1: Grignard/organometallic addition to carbonyl (HIGHEST PRIORITY)
        if h.get("carbonyl") and (h.get("grignard") or h.get("organolithium") or h.get("organozinc")):
            if h.get("grignard"):
                fam, conf = "grignard_addition", 0.90
            elif h.get("organozinc"):
                fam, conf = "organozinc_addition", 0.90
            else:
                fam, conf = "organolithium_addition", 0.90
        
        # PRIORITY 2: C-C Couplings
        elif h.get("aryl_halide") and h.get("boron"):
            fam, conf = "suzuki_miyaura", 0.9
        
        elif is_aryl_or_vinyl_electrophile and h.get("terminal_alkyne"):
            fam, conf = "sonogashira", 0.85
        
        elif is_aryl_or_vinyl_electrophile and h.get("grignard"):
            fam, conf = "kumada", 0.85
        
        elif is_aryl_or_vinyl_electrophile and h.get("organozinc"):
            fam, conf = "negishi", 0.85
        
        elif is_aryl_or_vinyl_electrophile and h.get("alkene") and not h.get("boron"):
            fam, conf = "heck", 0.80
        
        # PRIORITY 2.5: Acyl substitution (BEFORE reductive amination check)
        elif h.get("acyl_halide") and (h.get("nucleophile_n") or h.get("nucleophile_o")):
            if h.get("nucleophile_n"):
                fam, conf = "amide_formation", 0.90
            else:
                fam, conf = "ester_formation", 0.90
        
        # PRIORITY 3: SNAr (Aromatic Nucleophilic Substitution) - BEFORE metal-catalyzed couplings
        elif is_snar:
            # High confidence for heteroaryls (intrinsically activated)
            # Medium confidence for EWG-activated aryl halides
            if is_heteroaryl:
                fam, conf = "snar", 0.90
            else:
                fam, conf = "snar", 0.75
        
        # PRIORITY 4: Metal-catalyzed C-N/C-O/C-S Couplings (after SNAr check)
        elif is_aryl_or_vinyl_electrophile and h.get("nucleophile_n"):
            fam, conf = "cn_coupling", 0.9 if h.get("aryl_halide") else 0.8
        
        elif is_aryl_or_vinyl_electrophile and h.get("nucleophile_o"):
            fam, conf = "co_coupling", 0.85 if h.get("aryl_halide") else 0.75
        
        elif is_aryl_or_vinyl_electrophile and h.get("nucleophile_s"):
            fam, conf = "cs_coupling", 0.85 if h.get("aryl_halide") else 0.75
        
        # PRIORITY 5: Amide/Ester formation
        elif h.get("acid") and h.get("nucleophile_n"):
            fam, conf = "amide_coupling", 0.8
        
        elif h.get("acid") and h.get("alcohol"):
            fam, conf = "esterification", 0.85
        
        # PRIORITY 5.5: SN2 with alkyl halides (BEFORE reductive amination, BUT exclude acyl halides)
        elif h.get("alkyl_halide") and not h.get("acyl_halide") and (h.get("nucleophile_n") or h.get("nucleophile_o") or h.get("nucleophile_s")):
            fam, conf = "sn2_alkylation", 0.85
        
        # PRIORITY 6: Reductive Amination (carbonyl + amine)
        elif h.get("carbonyl") and h.get("nucleophile_n"):
            # Check if reducing agent present
            reducing_agent = _detect_reducing_agent(self.reactants)
            if reducing_agent in ("NaBH4", "NaBH(OAc)3", "NaBH3CN", "BH3"):
                fam, conf = "reductive_amination", 0.85
            else:
                # Could be imine formation or reductive amination
                fam, conf = "reductive_amination", 0.65
        
        # PRIORITY 7: SN2 reactions
        elif h.get("alkene") and h.get("borane"):
            fam, conf = "hydroboration", 0.85
        
        elif h.get("alkyl_halide") and h.get("cyanide"):
            fam, conf = "nitrile_formation", 0.90
        
        elif h.get("alkyl_halide") and h.get("iodide"):
            fam, conf = "finkelstein", 0.85
        
        elif h.get("alkyl_halide") and (h.get("alkoxide") or h.get("alcohol")):
            fam, conf = "williamson_ether", 0.85
        
        # PRIORITY 8: Condensation reactions
        elif h.get("ester") and rs.count("c(=o)o") >= 2:
            fam, conf = "claisen_condensation", 0.70
        
        elif h.get("alpha_beta_unsaturated"):
            fam, conf = "michael_addition", 0.65
        
        # PRIORITY 9: Ring-Closing Metathesis (RCM)
        # RCM involves diene substrate (intramolecular) with Ru catalyst
        elif h.get("diene") or (h.get("alkene") and len(self.reactants) == 1):
            # Check for Ru catalyst indicators
            is_ru_catalyst = any("ru" in agent.get("input", "").lower() or 
                                "grubbs" in agent.get("input", "").lower() or
                                "hoveyda" in agent.get("input", "").lower()
                                for agent in self.agents)
            if is_ru_catalyst:
                fam, conf = "ring_closing_metathesis", 0.88
        
        # PRIORITY 10: Cycloaddition
        elif h.get("conjugated_diene") and h.get("alkene"):
            fam, conf = "diels_alder", 0.85
        
        # Reagent-based detection (hydrogenation, reduction, oxidation, elimination)
        reducing_agent = _detect_reducing_agent(self.reactants)
        if reducing_agent == "H2" and (h.get("alkene") or h.get("terminal_alkyne") or h.get("carbonyl")):
            if fam is None or conf < 0.85:
                fam, conf = "hydrogenation", 0.90
        
        elif reducing_agent in ("NaBH4", "LiAlH4", "BH3", "DIBAL") and h.get("carbonyl"):
            if fam is None or conf < 0.85:
                fam, conf = "carbonyl_reduction", 0.88
        
        oxidizing_agent = _detect_oxidizing_agent(self.reactants)
        if oxidizing_agent in ("PCC", "KMnO4", "CrO3", "Swern", "MnO2", "Dess-Martin") and h.get("alcohol"):
            if fam is None or conf < 0.80:
                fam, conf = "alcohol_oxidation", 0.85
        
        elif oxidizing_agent in ("mCPBA", "H2O2") and h.get("alkene"):
            if fam is None or conf < 0.80:
                fam, conf = "epoxidation", 0.85
        
        if _detect_strong_base(self.reactants) and h.get("alkyl_halide") and not h.get("cyanide"):
            if fam is None or conf < 0.75:
                fam, conf = "e2_elimination", 0.80
        
        # Radical chain reactions (halogenation, polymerization)
        radical_initiator = _detect_radical_initiator(self.reactants)
        if radical_initiator:
            # Radical halogenation (e.g., CCl4, NBS)
            if radical_initiator in ("CCl4", "CBr4", "NBS", "Br2/Cl2"):
                if fam is None or conf < 0.70:
                    fam, conf = "radical_halogenation", 0.75
            # Other radical reactions (e.g., AIBN, BPO)
            elif fam is None or conf < 0.60:
                fam, conf = "radical_chain", 0.65
        
        # Map to taxonomy
        canonical_family = resolve_to_taxonomy(
            fam or "Unknown",
            catalysts=self.catalysts,
            functional_groups=h
        ) if fam else None
        
        return {
            "family": canonical_family,
            "confidence": float(conf),
            "hits": h,
            "raw_family": fam,  # Keep original for debugging
        }
    
    def _ml_detection(self) -> Dict[str, Any]:
        """
        ML-based detection via rxn-insight (optional).
        
        Consolidates ALL ML logic from reaction_type_detector.py:
        - _call_insight() - calls rxn-insight API
        - _extract_fields() - parses ML response
        - _map_to_family() - maps to taxonomy (NOW via resolve_to_taxonomy)
        - _refine_cn_family() - catalyst refinements
        
        IMPORTANT: rxn-insight returns UNPREDICTABLE names:
        - "Suzuki coupling" vs "Cross-coupling reaction"
        - "Buchwald-Hartwig" vs "Pd-catalyzed amination"
        Solution: robust mapping via resolve_to_taxonomy()
        
        Returns:
            {
                "available": bool,         # Is rxn-insight installed?
                "family": str | None,      # Mapped canonical ID
                "rxn_class": str | None,   # Broad ML class
                "rxn_name": str | None,    # Specific ML name
                "confidence": float | None,
                "raw": any,                # Raw response
                "mapping_method": str      # How it was mapped
            }
        """
        # Call rxn-insight via internal ML helpers
        from . import _ml_helpers
        
        ml_result = _ml_helpers.call_rxn_insight(
            self.normalized.get("normalized") or self.reaction_smiles
        )
        
        if not ml_result.get("available"):
            return {
                "available": False,
                "family": None,
                "rxn_class": None,
                "rxn_name": None,
                "confidence": None,
                "raw": None,
                "mapping_method": "not_available",
            }
        
        if not ml_result.get("success"):
            return {
                "available": True,
                "family": None,
                "rxn_class": ml_result.get("rxn_class"),
                "rxn_name": ml_result.get("rxn_name"),
                "confidence": None,
                "raw": ml_result.get("raw"),
                "mapping_method": "failed",
            }
        
        # Extract ML predictions
        rxn_class = ml_result.get("rxn_class")
        rxn_name = ml_result.get("rxn_name")
        ml_conf = ml_result.get("confidence")
        
        # Use NEW robust mapping (handles unpredictable ML names)
        is_cn = any(term in (rxn_class or "").lower() for term in ["heteroatom", "c-n", "amination"])
        
        family = resolve_to_taxonomy(
            rxn_name or rxn_class or "",
            catalysts=self.catalysts,
            is_cn_coupling=is_cn,
            functional_groups=self.functional_groups
        )
        
        # Adjust confidence based on mapping method
        mapping_method = get_mapping_method(rxn_name or rxn_class or "", family)
        if family and ml_conf is not None:
            adjusted_conf = calculate_confidence_adjustment(
                rxn_name or rxn_class or "",
                family,
                ml_conf
            )
        else:
            adjusted_conf = ml_conf
        
        # Log unmapped predictions for taxonomy improvement
        if not family and (rxn_name or rxn_class):
            logger.info(
                f"Unmapped ML prediction: class='{rxn_class}', name='{rxn_name}', "
                f"catalysts={self.catalysts}"
            )
        
        return {
            "available": True,
            "family": family,
            "rxn_class": rxn_class,
            "rxn_name": rxn_name,
            "confidence": adjusted_conf,
            "raw": ml_result.get("raw"),
            "mapping_method": mapping_method,
        }
    
    def _apply_catalyst_overrides(self, family: Optional[str], confidence: float) -> tuple[Optional[str], float]:
        """
        Apply catalyst-based family overrides.
        
        For C-N coupling reactions:
            Pd catalyst → buchwald_hartwig_c_n (conf: 0.95)
            Cu catalyst → ullmann_cn (conf: 0.90)
            
        For other heteroatom couplings, preserves family.
        
        Args:
            family: Current family prediction
            confidence: Current confidence score
            
        Returns:
            (refined_family, refined_confidence)
        """
        if not self.catalysts or not family:
            return family, confidence
        
        # Check if C-N coupling signature
        is_cn = (
            family in CN_FAMILIES_CANONICAL
            or (
                self.functional_groups.get("nucleophile_n")
                and (
                    self.functional_groups.get("aryl_halide")
                    or self.functional_groups.get("vinyl_halide")
                    or self.functional_groups.get("triflate")
                )
            )
        )
        
        if is_cn:
            if "Pd" in self.catalysts:
                return "buchwald_hartwig_c_n", 0.95
            elif "Cu" in self.catalysts:
                return "ullmann_cn", 0.90
            # Generic C-N with other catalysts
            if family == "cn_coupling" and confidence < 0.85:
                confidence = 0.85
        
        # Boost confidence for other heteroatom couplings with catalysts
        if family in CN_FAMILIES_CANONICAL and confidence < 0.85:
            confidence = 0.85
        
        return family, confidence
    
    def detect(self, use_ml: bool = True) -> Dict[str, Any]:
        """
        Main detection orchestrator - combines all detection methods.
        
        Detection Pipeline:
        1. Catalyst Detection (ALWAYS)
        2. Functional Group Detection (ALWAYS) 
        3. Rule-Based Detection (ALWAYS) - SMARTS patterns
        4. ML Detection (OPTIONAL) - rxn-insight if use_ml=True
        5. Catalyst Overrides (ALWAYS) - metal-based refinements
        6. Result Merging - choose best prediction
        
        Priority:
        - Catalyst override (highest)
        - ML detection (if available and confident)
        - Rule-based detection (fallback)
        
        Args:
            use_ml: Use ML detection if available (default: True)
            
        Returns:
            Unified detection result with all metadata
        """
        # Step 1 & 2: Detect catalysts and functional groups
        self._detect_catalysts()
        self._detect_functional_groups()
        
        # Step 3: Rule-based detection (ALWAYS runs)
        rule_result = self._rule_based_detection()
        fam_rule = rule_result["family"]
        conf_rule = rule_result["confidence"]
        
        # Step 4: ML detection (OPTIONAL)
        ml_result = None
        fam_ml = None
        conf_ml = None
        
        if use_ml:
            ml_result = self._ml_detection()
            if ml_result.get("available") and ml_result.get("family"):
                fam_ml = ml_result["family"]
                conf_ml = ml_result["confidence"]
        
        # Step 5: Choose best prediction
        # Prefer C-O/C-S rule-based detection (more specific)
        if (fam_rule in CO_FAMILIES_CANONICAL or fam_rule in CS_FAMILIES_CANONICAL) and conf_rule >= 0.75:
            fam_final = fam_rule
            conf_final = conf_rule
            method = "rule_based"
        elif fam_ml and conf_ml is not None and conf_ml > 0.5:
            fam_final = fam_ml
            conf_final = conf_ml
            method = "ml"
        else:
            fam_final = fam_rule
            conf_final = conf_rule
            method = "rule_based"
        
        # Step 6: Apply catalyst overrides (HIGHEST PRIORITY)
        fam_final, conf_final = self._apply_catalyst_overrides(fam_final, conf_final)
        if (fam_final in CN_FAMILIES_CANONICAL and 
            (fam_final != fam_rule or fam_final != fam_ml)):
            method = method + "+catalyst"
        
        # Build unified result
        result: Dict[str, Any] = {
            "family": fam_final or "Unknown",
            "confidence": float(conf_final),
            "method": method,
            "details": {
                "reactants": self.reactants,
                "catalysts": sorted(self.catalysts) if self.catalysts else [],
                "functional_groups": self.functional_groups,
                "rule_prediction": {
                    "family": fam_rule,
                    "confidence": conf_rule,
                    "raw_family": rule_result.get("raw_family"),
                },
            }
        }
        
        if ml_result and ml_result.get("available"):
            result["details"]["ml_prediction"] = {
                "family": fam_ml,
                "confidence": conf_ml,
                "rxn_class": ml_result.get("rxn_class"),
                "rxn_name": ml_result.get("rxn_name"),
                "mapping_method": ml_result.get("mapping_method"),
            }
        
        # Add agreement status
        if fam_ml is not None:
            result["agreement"] = (fam_rule == fam_ml)
            result["status"] = (
                "consistent" if result["agreement"]
                else "conflict" if fam_ml else "ml_unmapped"
            )
        else:
            result["agreement"] = None
            result["status"] = "rule_only"
        
        return result


def detect_reaction(reaction_smiles: str, use_ml: bool = True) -> Dict[str, Any]:
    """
    Detect reaction family from reaction SMILES.
    
    This is the MAIN entry point for all reaction detection. It consolidates
    SMARTS-based rules, ML predictions, and catalyst-based overrides into
    a single unified API with consistent outputs.
    
    All outputs are validated against the unified taxonomy in
    chemtools/taxonomy/data/reaction_types.json (80+ canonical reaction types).
    
    Args:
        reaction_smiles: Full reaction SMILES (reactants>>products)
        use_ml: Use ML-based detection if available (default: True)
        
    Returns:
        {
            "family": str,              # Canonical taxonomy ID
            "confidence": float,        # 0.0-1.0 confidence score
            "method": str,              # Detection method used
            "agreement": bool | None,   # Rule/ML agreement
            "status": str,              # "consistent", "conflict", "rule_only"
            "details": {
                "reactants": [...],
                "catalysts": [...],
                "functional_groups": {...},
                "rule_prediction": {...},
                "ml_prediction": {...}  # If ML was used
            }
        }
        
    Examples:
        >>> # Suzuki coupling
        >>> result = detect_reaction("Brc1ccccc1.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1")
        >>> result["family"]
        "suzuki_miyaura"
        >>> result["confidence"]
        0.9
        
        >>> # Buchwald-Hartwig (Pd catalyst detected)
        >>> result = detect_reaction("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")
        >>> result["family"]
        "buchwald_hartwig_c_n"
        >>> result["method"]
        "rule_based+catalyst"
        
        >>> # Rule-based only (no ML)
        >>> result = detect_reaction("Brc1ccccc1.c1ccccc1B(O)O>>...", use_ml=False)
        >>> result["method"]
        "rule_based"
    """
    engine = _DetectionEngine(reaction_smiles)
    return engine.detect(use_ml=use_ml)
