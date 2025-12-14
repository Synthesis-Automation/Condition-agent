"""
Unified reaction detection API.

This module consolidates all reaction detection logic (SMARTS-based, ML-based,
catalyst-based) into a single, clean API with consistent outputs.

Public Functions:
    detect_reaction(reaction_smiles, use_ml=True) -> dict
    
Internal:
    _DetectionEngine class (consolidates all detection logic)

Taxonomy-based detection (SMARTS matching, DRFP similarity) is handled by
the chemtools.taxonomy.detection module.

All outputs are validated against the unified taxonomy in chemtools/taxonomy/.
"""

from typing import Any, Dict, List, Optional, Set, Tuple
import logging
from pathlib import Path

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

# Import taxonomy-based detection functions from the dedicated module
from .taxonomy.detection import (
    detect_by_taxonomy_smarts,
    detect_by_drfp_similarity,
    _load_taxonomy_smarts,
    _load_reference_reactions,
    _compute_drfp_fingerprint,
    _get_reference_drfp_index,
    _validate_reaction_by_running,
    clear_caches as clear_taxonomy_caches,
)

logger = logging.getLogger(__name__)


# Note: Taxonomy-based detection functions (detect_by_taxonomy_smarts, detect_by_drfp_similarity,
# _load_taxonomy_smarts, _load_reference_reactions, _compute_drfp_fingerprint, _get_reference_drfp_index,
# _validate_reaction_by_running) are now imported from chemtools.taxonomy.detection module.


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
    reaction_smiles: str,
    threshold: float = 0.3,
    top_k: int = 3
) -> List[Tuple[str, str, float]]:
    """
    Detect reaction type by DRFP similarity to reference reactions.
    
    Computes the DRFP fingerprint of the query reaction and compares it
    to precomputed fingerprints of reference reactions in the taxonomy.
    
    Args:
        reaction_smiles: Query reaction SMILES
        threshold: Minimum similarity threshold (0.0 - 1.0)
        top_k: Maximum number of matches to return
        
    Returns:
        List of (reaction_type_id, name, similarity) tuples, sorted by similarity
    """
    try:
        import numpy as np
    except ImportError:
        logger.debug("NumPy not available for DRFP similarity")
        return []
    
    # Compute query DRFP
    query_drfp = _compute_drfp_fingerprint(reaction_smiles)
    if query_drfp is None:
        return []
    
    # Get reference DRFP index
    ref_index = _get_reference_drfp_index()
    if not ref_index:
        return []
    
    # Load taxonomy for names
    taxonomy = _load_taxonomy_smarts()
    
    # Compute similarities to all reference reactions
    matches = []
    
    for rxn_type, ref_list in ref_index.items():
        best_sim = 0.0
        for ref_smiles, ref_drfp in ref_list:
            # Cosine similarity
            dot = np.dot(query_drfp, ref_drfp)
            norm_q = np.linalg.norm(query_drfp)
            norm_r = np.linalg.norm(ref_drfp)
            
            if norm_q > 0 and norm_r > 0:
                sim = dot / (norm_q * norm_r)
                best_sim = max(best_sim, sim)
        
        if best_sim >= threshold:
            name = taxonomy.get(rxn_type, {}).get("name", rxn_type)
            matches.append((rxn_type, name, best_sim))
    
    # Sort by similarity (highest first)
    matches.sort(key=lambda x: x[2], reverse=True)
    
    return matches[:top_k]


def _load_taxonomy_smarts() -> Dict[str, Any]:
    """
    Load SMARTS patterns from reaction_types.json taxonomy file.
    
    Returns:
        Dict mapping reaction_type_id to {smarts, name, category}
    """
    global _TAXONOMY_SMARTS_CACHE
    
    if _TAXONOMY_SMARTS_CACHE is not None:
        return _TAXONOMY_SMARTS_CACHE
    
    import json
    
    taxonomy_path = Path(__file__).parent / "taxonomy" / "data" / "reaction_types.json"
    
    if not taxonomy_path.exists():
        logger.warning(f"Taxonomy file not found: {taxonomy_path}")
        _TAXONOMY_SMARTS_CACHE = {}
        return _TAXONOMY_SMARTS_CACHE
    
    try:
        with open(taxonomy_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        _TAXONOMY_SMARTS_CACHE = {}
        for entry in data:
            rxn_id = entry.get("id")
            smarts = entry.get("smarts")
            if rxn_id and smarts:
                _TAXONOMY_SMARTS_CACHE[rxn_id] = {
                    "smarts": smarts,
                    "name": entry.get("name", rxn_id),
                    "category": entry.get("category", ""),
                    "aliases": entry.get("aliases", []),
                }
        
        logger.debug(f"Loaded {len(_TAXONOMY_SMARTS_CACHE)} reaction SMARTS from taxonomy")
        return _TAXONOMY_SMARTS_CACHE
        
    except Exception as e:
        logger.warning(f"Failed to load taxonomy SMARTS: {e}")
        _TAXONOMY_SMARTS_CACHE = {}
        return _TAXONOMY_SMARTS_CACHE


def _validate_reaction_by_running(
    rxn_pattern,
    reactant_mols: List,
    product_mols: List,
    max_product_sets: int = 50
) -> Tuple[bool, float]:
    """
    Validate a reaction by running the SMARTS and comparing predicted vs actual products.
    
    This is more rigorous than just checking substructure matches because it:
    1. Runs the reaction SMARTS on reactant molecules
    2. Generates predicted product structures
    3. Compares predicted products with actual products
    
    Args:
        rxn_pattern: RDKit ChemicalReaction object
        reactant_mols: List of reactant RDKit Mol objects
        product_mols: List of actual product RDKit Mol objects
        max_product_sets: Maximum number of product sets to check
        
    Returns:
        (is_match, confidence_boost)
        - is_match: True if predicted products match actual products
        - confidence_boost: Additional confidence based on match quality
    """
    try:
        from rdkit import Chem
        from itertools import permutations
        
        if not reactant_mols or not product_mols:
            return False, 0.0
        
        num_templates = rxn_pattern.GetNumReactantTemplates()
        if num_templates == 0:
            return False, 0.0
        
        # Get canonical SMILES of actual products for comparison
        actual_product_smiles = set()
        for mol in product_mols:
            try:
                smi = Chem.MolToSmiles(mol, canonical=True)
                actual_product_smiles.add(smi)
            except Exception:
                continue
        
        if not actual_product_smiles:
            return False, 0.0
        
        # Try different orderings of reactants (in case order matters)
        # Limit to reasonable number of permutations
        reactant_orderings = []
        if len(reactant_mols) <= 4:
            for perm in permutations(range(len(reactant_mols))):
                if len(perm) >= num_templates:
                    reactant_orderings.append([reactant_mols[i] for i in perm[:num_templates]])
                if len(reactant_orderings) >= 10:  # Limit permutations
                    break
        else:
            # Just use first N reactants
            reactant_orderings.append(reactant_mols[:num_templates])
        
        # Run reaction and check if any predicted products match actual products
        best_match_score = 0.0
        
        for reactants in reactant_orderings:
            if len(reactants) != num_templates:
                continue
            
            try:
                product_sets = rxn_pattern.RunReactants(tuple(reactants))
            except Exception:
                continue
            
            # Check each product set
            for prod_set in product_sets[:max_product_sets]:
                predicted_smiles = set()
                for pred_mol in prod_set:
                    try:
                        # Sanitize the predicted product
                        Chem.SanitizeMol(pred_mol)
                        smi = Chem.MolToSmiles(pred_mol, canonical=True)
                        predicted_smiles.add(smi)
                    except Exception:
                        continue
                
                if not predicted_smiles:
                    continue
                
                # Calculate overlap between predicted and actual products
                # Match score = intersection / union (Jaccard similarity)
                intersection = predicted_smiles & actual_product_smiles
                union = predicted_smiles | actual_product_smiles
                
                if intersection:
                    match_score = len(intersection) / len(union)
                    best_match_score = max(best_match_score, match_score)
                    
                    # If we found a perfect or near-perfect match, return early
                    if match_score >= 0.8:
                        return True, match_score * 0.15  # Up to +0.15 confidence boost
        
        # Return result based on best match
        if best_match_score >= 0.5:
            return True, best_match_score * 0.15
        
        return False, 0.0
        
    except Exception as e:
        logger.debug(f"Reaction validation failed: {e}")
        return False, 0.0


def detect_by_taxonomy_smarts(reaction_smiles: str) -> List[Tuple[str, str, float]]:
    """
    Detect reaction type by matching against taxonomy SMARTS patterns.
    
    This function uses the SMARTS patterns defined in reaction_types.json
    to identify the reaction type. Each pattern is a reaction SMARTS that
    describes the transformation.
    
    Args:
        reaction_smiles: Full reaction SMILES (reactants>>products)
        
    Returns:
        List of (reaction_type_id, name, confidence) tuples, sorted by confidence.
        Empty list if no matches found.
        
    Example:
        >>> matches = detect_by_taxonomy_smarts("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
        >>> matches[0]
        ('suzuki_miyaura', 'Suzuki-Miyaura Coupling', 0.95)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.warning("RDKit not available for taxonomy SMARTS matching")
        return []
    
    taxonomy = _load_taxonomy_smarts()
    if not taxonomy:
        return []
    
    matches: List[Tuple[str, str, float]] = []
    
    # Parse reaction SMILES
    try:
        rxn_mol = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        if rxn_mol is None:
            # Try parsing as regular reaction
            parts = reaction_smiles.split(">>")
            if len(parts) != 2:
                return []
            reactants_smiles, products_smiles = parts
    except Exception:
        # Fallback: parse reactants and products separately
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            return []
        reactants_smiles, products_smiles = parts
    
    # Parse reactants and products as molecules
    try:
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            return []
        
        reactants_smiles, products_smiles = parts
        
        # Parse individual reactants
        reactant_mols = []
        for smi in reactants_smiles.split("."):
            mol = Chem.MolFromSmiles(smi)
            if mol:
                reactant_mols.append(mol)
        
        # Parse products
        product_mols = []
        for smi in products_smiles.split("."):
            mol = Chem.MolFromSmiles(smi)
            if mol:
                product_mols.append(mol)
        
        if not reactant_mols or not product_mols:
            return []
            
    except Exception as e:
        logger.debug(f"Failed to parse reaction: {e}")
        return []
    
    # Try matching each taxonomy SMARTS
    for rxn_id, info in taxonomy.items():
        smarts = info["smarts"]
        
        try:
            # Compile reaction SMARTS
            rxn_pattern = AllChem.ReactionFromSmarts(smarts)
            if rxn_pattern is None:
                continue
            
            # Check if reactants AND products match the pattern
            try:
                # Get number of reactant and product templates
                num_reactant_templates = rxn_pattern.GetNumReactantTemplates()
                num_product_templates = rxn_pattern.GetNumProductTemplates()
                
                if num_reactant_templates == 0:
                    continue
                
                # === STEP 1: Match reactant templates ===
                # Each template must match at least one reactant
                template_matches = []  # List of sets: which reactants each template matches
                
                for i in range(num_reactant_templates):
                    template = rxn_pattern.GetReactantTemplate(i)
                    matching_reactants = set()
                    for j, mol in enumerate(reactant_mols):
                        if mol.HasSubstructMatch(template):
                            matching_reactants.add(j)
                    template_matches.append(matching_reactants)
                
                # All reactant templates must have at least one match
                if any(len(m) == 0 for m in template_matches):
                    continue
                
                # Try to find a valid assignment (simple greedy for 2 templates)
                reactants_matched = False
                if num_reactant_templates == 1:
                    reactants_matched = len(template_matches[0]) > 0
                elif num_reactant_templates == 2:
                    # Check if we can assign different reactants to each template
                    for r1 in template_matches[0]:
                        for r2 in template_matches[1]:
                            if r1 != r2:
                                reactants_matched = True
                                break
                        if reactants_matched:
                            break
                    # Also accept if we have exactly 2 reactants and each matches its template
                    if not reactants_matched and len(reactant_mols) == 2:
                        if 0 in template_matches[0] and 1 in template_matches[1]:
                            reactants_matched = True
                        elif 1 in template_matches[0] and 0 in template_matches[1]:
                            reactants_matched = True
                else:
                    # For 3+ templates, use simple check: all have matches
                    reactants_matched = all(len(m) > 0 for m in template_matches)
                
                if not reactants_matched:
                    continue
                
                # === STEP 2: Match product templates (NEW) ===
                # Each product template must match at least one product
                products_matched = True  # Default to True if no product templates
                
                if num_product_templates > 0 and product_mols:
                    product_template_matches = []
                    
                    for i in range(num_product_templates):
                        template = rxn_pattern.GetProductTemplate(i)
                        matching_products = set()
                        for j, mol in enumerate(product_mols):
                            if mol.HasSubstructMatch(template):
                                matching_products.add(j)
                        product_template_matches.append(matching_products)
                    
                    # All product templates must have at least one match
                    if any(len(m) == 0 for m in product_template_matches):
                        products_matched = False
                    else:
                        # For multiple product templates, check valid assignment
                        if num_product_templates == 1:
                            products_matched = len(product_template_matches[0]) > 0
                        elif num_product_templates == 2:
                            # Check if we can assign different products to each template
                            products_matched = False
                            for p1 in product_template_matches[0]:
                                for p2 in product_template_matches[1]:
                                    if p1 != p2:
                                        products_matched = True
                                        break
                                if products_matched:
                                    break
                            # Also accept if exactly 2 products match their templates
                            if not products_matched and len(product_mols) == 2:
                                if 0 in product_template_matches[0] and 1 in product_template_matches[1]:
                                    products_matched = True
                                elif 1 in product_template_matches[0] and 0 in product_template_matches[1]:
                                    products_matched = True
                        else:
                            products_matched = all(len(m) > 0 for m in product_template_matches)
                
                # === STEP 3: Only match if BOTH reactants AND products match ===
                if reactants_matched and products_matched:
                    # Calculate confidence based on specificity
                    # More specific patterns (longer SMARTS) get higher confidence
                    specificity = min(len(smarts) / 150.0, 1.0)  # Adjusted scale
                    
                    # Bonus for exact number of reactant match
                    if len(reactant_mols) == num_reactant_templates:
                        specificity = min(specificity + 0.1, 1.0)
                    
                    # Bonus for product validation (more reliable match)
                    if num_product_templates > 0 and products_matched:
                        specificity = min(specificity + 0.05, 1.0)
                    
                    # === STEP 4: Run reaction to validate predicted products match actual ===
                    # This is the most rigorous validation - actually run the reaction
                    reaction_validated, validation_boost = _validate_reaction_by_running(
                        rxn_pattern, reactant_mols, product_mols
                    )
                    
                    if reaction_validated:
                        # Strong boost for reactions that produce matching products
                        specificity = min(specificity + validation_boost, 1.0)
                    
                    confidence = 0.75 + (0.20 * specificity)
                    
                    matches.append((rxn_id, info["name"], confidence, reaction_validated))
                    
            except Exception:
                # Reaction running failed, try simple substructure matching
                continue
                
        except Exception as e:
            logger.debug(f"Failed to match SMARTS for {rxn_id}: {e}")
            continue
            continue
    
    # Sort by: 1) reaction_validated (True first), 2) confidence (highest first)
    matches.sort(key=lambda x: (x[3] if len(x) > 3 else False, x[2]), reverse=True)
    
    # Return without the validation flag (keep API compatible)
    return [(m[0], m[1], m[2]) for m in matches]


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
        self.products: List[str] = []  # NEW: Track products for validation
        self.agents: List[Dict[str, Any]] = []
        self.catalysts: Set[str] = set()
        self.functional_groups: Dict[str, bool] = {}
        self.product_functional_groups: Dict[str, bool] = {}  # NEW: Product FGs
        
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
            
            # NEW: Extract products for validation
            self.products = [
                (p.get("smiles_norm") or p.get("largest_smiles") or p.get("input") or "")
                for p in (self.normalized.get("products") or [])
            ]
            self.products = [s for s in self.products if s]
            
            # Extract agents
            self.agents = self.normalized.get("agents") or []
            
        except Exception as e:
            logger.warning(f"Failed to normalize reaction: {e}")
            self.normalized = {}
            self.reactants = []
            self.products = []
            self.agents = []
    
    def _detect_catalysts(self) -> Set[str]:
        """
        Detect catalyst metals from agents AND reactants.
        
        Extracts Pd, Cu, Ni, Co, Ru, Rh, Ir, Fe from:
        - Agent position (between > and >)
        - Reactant position (some reactions have catalyst in reactants)
        
        Returns:
            Set of metal symbols (e.g., {"Pd", "Cu"})
        """
        try:
            self.catalysts = _detect_agent_metals(self.agents)
        except Exception as e:
            logger.warning(f"Failed to detect catalysts from agents: {e}")
            self.catalysts = set()
        
        # NEW: Also check reactants for metal catalysts
        # (Some reactions encode catalyst in reactant position)
        metal_patterns = {
            "Pd": ["[Pd]", "Pd(", "pd(", "palladium"],
            "Cu": ["[Cu]", "Cu(", "cu(", "copper", "CuI", "CuBr", "CuCl"],
            "Ni": ["[Ni]", "Ni(", "ni(", "nickel"],
            "Co": ["[Co]", "Co(", "co(", "cobalt"],
            "Ru": ["[Ru]", "Ru(", "ru(", "ruthenium", "grubbs", "hoveyda"],
            "Rh": ["[Rh]", "Rh(", "rh(", "rhodium"],
            "Ir": ["[Ir]", "Ir(", "ir(", "iridium"],
            "Fe": ["[Fe]", "Fe(", "fe(", "iron", "FeCl"],
        }
        
        reactants_str = " ".join(self.reactants).lower()
        for metal, patterns in metal_patterns.items():
            if any(p.lower() in reactants_str for p in patterns):
                self.catalysts.add(metal)
        
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
        
        # NEW: Also detect functional groups in products for validation
        try:
            self.product_functional_groups = _rule_hits(self.products) if self.products else {}
        except Exception as e:
            logger.warning(f"Failed to detect product functional groups: {e}")
            self.product_functional_groups = {}
        
        # Add catalyst info to functional groups
        if self.catalysts:
            self.functional_groups["catalyst_pd"] = "Pd" in self.catalysts
            self.functional_groups["catalyst_cu"] = "Cu" in self.catalysts
            self.functional_groups["catalyst_ni"] = "Ni" in self.catalysts
            self.functional_groups["catalyst_co"] = "Co" in self.catalysts
            self.functional_groups["catalyst_ru"] = "Ru" in self.catalysts
            self.functional_groups["catalyst_rh"] = "Rh" in self.catalysts
        
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
        elif (h.get("aryl_halide") or h.get("vinyl_halide")) and h.get("boron"):
            fam, conf = "suzuki_miyaura", 0.9
        
        elif is_aryl_or_vinyl_electrophile and h.get("terminal_alkyne"):
            fam, conf = "sonogashira", 0.85
        
        elif is_aryl_or_vinyl_electrophile and h.get("grignard"):
            fam, conf = "kumada", 0.85
        
        elif is_aryl_or_vinyl_electrophile and h.get("organozinc"):
            fam, conf = "negishi", 0.85
        
        elif is_aryl_or_vinyl_electrophile and h.get("organostannane"):
            fam, conf = "stille", 0.85
        
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
    
    def _taxonomy_smarts_detection(self) -> Dict[str, Any]:
        """
        Taxonomy SMARTS-based detection using reaction_types.json SMARTS patterns.
        
        Uses the SMARTS patterns defined in reaction_types.json to match the reaction
        against known reaction type signatures.
        
        Returns:
            {
                "matches": list,           # List of (reaction_type_id, name, confidence)
                "best_family": str | None, # Best matching reaction type
                "confidence": float | None,
            }
        """
        try:
            matches = detect_by_taxonomy_smarts(self.reaction_smiles)
            if matches:
                best = matches[0]  # Already sorted by confidence
                return {
                    "matches": matches,
                    "best_family": best[0],  # reaction_type_id
                    "confidence": best[2],   # confidence score
                }
        except Exception as e:
            logger.warning(f"Taxonomy SMARTS detection failed: {e}")
        
        return {
            "matches": [],
            "best_family": None,
            "confidence": None,
        }
    
    def _drfp_similarity_detection(self) -> Dict[str, Any]:
        """
        DRFP similarity-based detection using reference reactions.
        
        Computes DRFP fingerprint of the query reaction and compares to
        precomputed fingerprints of reference reactions in the taxonomy.
        This provides a complementary ML-based detection method.
        
        Returns:
            {
                "matches": list,           # List of (reaction_type_id, name, similarity)
                "best_family": str | None, # Best matching reaction type
                "confidence": float | None,
                "available": bool,         # Whether DRFP is available
            }
        """
        try:
            matches = detect_by_drfp_similarity(
                self.reaction_smiles,
                threshold=0.3,
                top_k=5
            )
            if matches:
                best = matches[0]  # Already sorted by similarity
                # Convert similarity (0-1) to confidence
                # DRFP similarity of 0.5+ is quite good for reaction matching
                confidence = min(best[2] * 1.2, 1.0)  # Slight boost, cap at 1.0
                return {
                    "matches": matches,
                    "best_family": best[0],
                    "confidence": confidence,
                    "available": True,
                }
        except Exception as e:
            logger.debug(f"DRFP similarity detection failed: {e}")
        
        return {
            "matches": [],
            "best_family": None,
            "confidence": None,
            "available": False,
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
    
    def _validate_with_product(self, family: str, confidence: float) -> Tuple[str, float]:
        """
        NEW: Validate/adjust prediction based on product analysis.
        
        Checks if product structure is consistent with predicted reaction type:
        - C-N coupling → product should have C-N bond (no N nucleophile remaining)
        - Suzuki → product should be biaryl (no boron remaining)
        - Reduction → carbonyl consumed
        
        Args:
            family: Predicted family
            confidence: Current confidence
            
        Returns:
            (validated_family, adjusted_confidence)
        """
        if not self.product_functional_groups or not family:
            return family, confidence
        
        r_fg = self.functional_groups
        p_fg = self.product_functional_groups
        
        # C-N Coupling validation: 
        # The nitrogen is NOT consumed - it's still in the product, bonded to aryl carbon
        # We check: 1) aryl halide consumed, 2) product has aryl-N bond
        if family in CN_FAMILIES_CANONICAL or family == "cn_coupling":
            # Aryl halide should be consumed
            halide_consumed = r_fg.get("aryl_halide") and not p_fg.get("aryl_halide")
            if halide_consumed:
                confidence = min(confidence + 0.05, 1.0)
            # Note: We can't easily check for aryl-N bond formation without more complex analysis
        
        # Suzuki validation: boron should be consumed (not present in product)
        if family == "suzuki_miyaura":
            if r_fg.get("boron") and not p_fg.get("boron"):
                # Boron consumed - confirms Suzuki
                confidence = min(confidence + 0.05, 1.0)
            elif r_fg.get("boron") and p_fg.get("boron"):
                # Boron still present - lower confidence
                confidence = max(confidence - 0.1, 0.5)
        
        # Aryl/vinyl halide should be consumed in coupling reactions
        coupling_families = {"suzuki_miyaura", "sonogashira", "heck", "negishi", 
                            "kumada", "stille", "buchwald_hartwig_c_n", "ullmann_cn",
                            "cn_coupling", "co_coupling", "cs_coupling"}
        if family in coupling_families:
            halide_consumed = (
                (r_fg.get("aryl_halide") and not p_fg.get("aryl_halide")) or
                (r_fg.get("vinyl_halide") and not p_fg.get("vinyl_halide"))
            )
            if halide_consumed:
                # Halide consumed - consistent with coupling
                confidence = min(confidence + 0.03, 1.0)
        
        # Hydrogenation validation: alkene should be consumed
        if family == "hydrogenation":
            if r_fg.get("alkene") and p_fg.get("alkene"):
                # Alkene still present in product - NOT hydrogenation
                # This is likely Heck or other reaction that preserves alkene
                confidence = max(confidence - 0.5, 0.0)  # Strong penalty
            elif r_fg.get("alkene") and not p_fg.get("alkene"):
                # Alkene consumed - confirms hydrogenation
                confidence = min(confidence + 0.1, 1.0)
        
        return family, confidence

    def detect(self, use_ml: bool = True, use_taxonomy_smarts: bool = True) -> Dict[str, Any]:
        """
        Main detection orchestrator - combines all detection methods.
        
        Detection Pipeline:
        1. Catalyst Detection (ALWAYS)
        2. Functional Group Detection (ALWAYS) 
        3. Rule-Based Detection (ALWAYS) - SMARTS patterns
        4. Taxonomy SMARTS Detection (OPTIONAL) - reaction_types.json patterns
        5. ML Detection (OPTIONAL) - rxn-insight if use_ml=True
        6. Catalyst Overrides (ALWAYS) - metal-based refinements
        7. Result Merging - choose best prediction
        
        Priority:
        - Catalyst override (highest)
        - Taxonomy SMARTS match (highest confidence, exact match)
        - ML detection (if available and confident)
        - Rule-based detection (fallback)
        
        Args:
            use_ml: Use ML detection if available (default: True)
            use_taxonomy_smarts: Use taxonomy SMARTS detection (default: True)
            
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
        
        # Step 3.5: Taxonomy SMARTS detection (OPTIONAL)
        taxonomy_result = None
        fam_taxonomy = None
        conf_taxonomy = None
        
        if use_taxonomy_smarts:
            taxonomy_result = self._taxonomy_smarts_detection()
            if taxonomy_result.get("best_family"):
                fam_taxonomy = taxonomy_result["best_family"]
                conf_taxonomy = taxonomy_result["confidence"]
        
        # Step 4: ML detection (OPTIONAL)
        ml_result = None
        fam_ml = None
        conf_ml = None
        
        if use_ml:
            ml_result = self._ml_detection()
            if ml_result.get("available") and ml_result.get("family"):
                fam_ml = ml_result["family"]
                conf_ml = ml_result["confidence"]
        
        # Step 4.5: DRFP similarity detection (NEW - uses reference reactions)
        # DRFP is fast (fingerprint comparison) and always runs
        drfp_result = None
        fam_drfp = None
        conf_drfp = None
        
        drfp_result = self._drfp_similarity_detection()
        if drfp_result.get("available") and drfp_result.get("best_family"):
            fam_drfp = drfp_result["best_family"]
            conf_drfp = drfp_result["confidence"]
        
        # Step 5: Choose best prediction
        # Priority: taxonomy_smarts (high conf) > C-O/C-S rule > ML > rule_based
        # EXCEPTION: If rule-based detects C-N/C-O/C-S coupling with nucleophile,
        # prefer it over taxonomy C-C coupling (e.g., avoid Negishi misclassification)
        
        # Check if rule-based detected a heteroatom coupling with clear evidence
        is_rule_heteroatom_coupling = (
            fam_rule in CN_FAMILIES_CANONICAL or 
            fam_rule in CO_FAMILIES_CANONICAL or 
            fam_rule in CS_FAMILIES_CANONICAL or
            fam_rule in ("cn_coupling", "co_coupling", "cs_coupling")
        ) and conf_rule >= 0.85
        
        # Check if taxonomy is suggesting a C-C coupling (potential misclassification)
        taxonomy_is_cc_coupling = (
            fam_taxonomy in ("negishi", "kumada", "suzuki_miyaura", "stille", "heck", "sonogashira")
        ) if fam_taxonomy else False
        
        # Check if DRFP detects a C-C coupling (useful for unusual substrates)
        drfp_is_cc_coupling = (
            fam_drfp in ("negishi", "kumada", "suzuki_miyaura", "stille", "heck", "sonogashira")
        ) if fam_drfp else False
        
        # Check if DRFP suggests Suzuki and boron is present (strong signal)
        drfp_suzuki_with_boron = (
            fam_drfp == "suzuki_miyaura" and 
            "boron" in self.functional_groups
        ) if fam_drfp else False
        
        # If rule says heteroatom coupling but taxonomy says C-C coupling, trust rule
        if is_rule_heteroatom_coupling and taxonomy_is_cc_coupling:
            fam_final = fam_rule
            conf_final = conf_rule
            method = "rule_based"
        # DRFP Suzuki with boron present beats generic taxonomy (e.g., sn2_substitution)
        # This handles unusual electrophiles like cyclopropyl bromide
        elif drfp_suzuki_with_boron and fam_drfp != fam_taxonomy and conf_drfp >= 0.4:
            fam_final = fam_drfp
            conf_final = conf_drfp + 0.1  # Boost confidence since boron confirms
            method = "drfp_similarity+boron"
        # Prefer C-C coupling rule-based detection (very specific patterns)
        # This protects against generic taxonomy patterns like hydrogenation/grignard
        elif fam_rule in ("suzuki_miyaura", "heck", "sonogashira", "negishi", "kumada", "stille") and conf_rule >= 0.75:
            fam_final = fam_rule
            conf_final = conf_rule
            method = "rule_based"
        # Then check taxonomy SMARTS (exact pattern match is very reliable)
        elif fam_taxonomy and conf_taxonomy is not None and conf_taxonomy >= 0.8:
            fam_final = fam_taxonomy
            conf_final = conf_taxonomy
            method = "taxonomy_smarts"
        # Prefer C-O/C-S rule-based detection (more specific)
        elif (fam_rule in CO_FAMILIES_CANONICAL or fam_rule in CS_FAMILIES_CANONICAL) and conf_rule >= 0.75:
            fam_final = fam_rule
            conf_final = conf_rule
            method = "rule_based"
        # Check ML detection (high confidence)
        elif fam_ml and conf_ml is not None and conf_ml > 0.7:
            fam_final = fam_ml
            conf_final = conf_ml
            method = "ml"
        # Check DRFP similarity detection (NEW - good for specific reaction types)
        elif fam_drfp and conf_drfp is not None and conf_drfp >= 0.7:
            fam_final = fam_drfp
            conf_final = conf_drfp
            method = "drfp_similarity"
        # Check ML detection with lower confidence
        elif fam_ml and conf_ml is not None and conf_ml > 0.5:
            fam_final = fam_ml
            conf_final = conf_ml
            method = "ml"
        # Check DRFP with moderate confidence
        elif fam_drfp and conf_drfp is not None and conf_drfp >= 0.5:
            fam_final = fam_drfp
            conf_final = conf_drfp
            method = "drfp_similarity"
        # Check taxonomy SMARTS with lower confidence
        elif fam_taxonomy and conf_taxonomy is not None and conf_taxonomy >= 0.5:
            fam_final = fam_taxonomy
            conf_final = conf_taxonomy
            method = "taxonomy_smarts"
        else:
            fam_final = fam_rule
            conf_final = conf_rule
            method = "rule_based"
        
        # Step 6: Apply catalyst overrides (HIGHEST PRIORITY)
        fam_final, conf_final = self._apply_catalyst_overrides(fam_final, conf_final)
        if (fam_final in CN_FAMILIES_CANONICAL and 
            (fam_final != fam_rule or fam_final != fam_ml)):
            method = method + "+catalyst"
        
        # Step 7: Validate with product analysis (NEW)
        if self.products:
            fam_final, conf_final = self._validate_with_product(fam_final, conf_final)
        
        # Build unified result
        result: Dict[str, Any] = {
            "family": fam_final or "Unknown",
            "confidence": float(conf_final),
            "method": method,
            "details": {
                "reactants": self.reactants,
                "products": self.products,  # NEW: Include products
                "catalysts": sorted(self.catalysts) if self.catalysts else [],
                "functional_groups": self.functional_groups,
                "product_functional_groups": self.product_functional_groups if self.product_functional_groups else None,  # NEW
                "rule_prediction": {
                    "family": fam_rule,
                    "confidence": conf_rule,
                    "raw_family": rule_result.get("raw_family"),
                },
            }
        }
        
        # Add taxonomy SMARTS prediction if available
        if taxonomy_result and taxonomy_result.get("best_family"):
            result["details"]["taxonomy_smarts_prediction"] = {
                "family": fam_taxonomy,
                "confidence": conf_taxonomy,
                "all_matches": taxonomy_result.get("matches", []),
            }
        
        if ml_result and ml_result.get("available"):
            result["details"]["ml_prediction"] = {
                "family": fam_ml,
                "confidence": conf_ml,
                "rxn_class": ml_result.get("rxn_class"),
                "rxn_name": ml_result.get("rxn_name"),
                "mapping_method": ml_result.get("mapping_method"),
            }
        
        # Add DRFP similarity prediction if available (NEW)
        if drfp_result and drfp_result.get("available"):
            result["details"]["drfp_prediction"] = {
                "family": fam_drfp,
                "confidence": conf_drfp,
                "matched_reference": drfp_result.get("matched_reference"),
                "similarity_score": drfp_result.get("similarity"),
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


def detect_reaction(reaction_smiles: str, use_ml: bool = True, use_taxonomy_smarts: bool = True) -> Dict[str, Any]:
    """
    Detect reaction family from reaction SMILES.
    
    This is the MAIN entry point for all reaction detection. It consolidates
    SMARTS-based rules, ML predictions, taxonomy SMARTS matching, and catalyst-based 
    overrides into a single unified API with consistent outputs.
    
    All outputs are validated against the unified taxonomy in
    chemtools/taxonomy/data/reaction_types.json (80+ canonical reaction types).
    
    Args:
        reaction_smiles: Full reaction SMILES (reactants>>products)
        use_ml: Use ML-based detection if available (default: True)
        use_taxonomy_smarts: Use taxonomy SMARTS matching (default: True)
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
                "taxonomy_smarts_prediction": {...},  # If taxonomy SMARTS matched
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
        
        >>> # Taxonomy SMARTS detection
        >>> result = detect_reaction("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
        >>> result["method"]
        "taxonomy_smarts"
    """
    engine = _DetectionEngine(reaction_smiles)
    return engine.detect(use_ml=use_ml, use_taxonomy_smarts=use_taxonomy_smarts)
