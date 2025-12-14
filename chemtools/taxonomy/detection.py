"""
Taxonomy-based reaction type detection.

This module provides reaction type detection using:
1. SMARTS pattern matching against taxonomy definitions
2. DRFP (Differential Reaction Fingerprint) similarity to reference reactions
3. Reaction running validation (comparing predicted vs actual products)

Public Functions:
    detect_by_taxonomy_smarts(reaction_smiles) -> List[Tuple[str, str, float]]
    detect_by_drfp_similarity(reaction_smiles, threshold, top_k) -> List[Tuple[str, str, float]]

The taxonomy is defined in chemtools/taxonomy/data/reaction_types.json.
"""

from typing import Any, Dict, List, Optional, Tuple
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ============================================================================
# Module-level caches (loaded once)
# ============================================================================
_TAXONOMY_SMARTS_CACHE: Optional[Dict[str, Any]] = None
_REFERENCE_REACTIONS_CACHE: Optional[Dict[str, List[str]]] = None
_REFERENCE_DRFP_CACHE: Optional[Dict[str, Any]] = None


def _get_taxonomy_path() -> Path:
    """Get path to reaction_types.json taxonomy file."""
    return Path(__file__).parent / "data" / "reaction_types.json"


# ============================================================================
# Reference Reactions Loading (for DRFP similarity)
# ============================================================================

def _load_reference_reactions() -> Dict[str, List[str]]:
    """
    Load reference reactions from reaction_types.json taxonomy file.
    
    Returns:
        Dict mapping reaction_type_id to list of reference reaction SMILES
    """
    global _REFERENCE_REACTIONS_CACHE
    
    if _REFERENCE_REACTIONS_CACHE is not None:
        return _REFERENCE_REACTIONS_CACHE
    
    import json
    
    taxonomy_path = _get_taxonomy_path()
    
    if not taxonomy_path.exists():
        _REFERENCE_REACTIONS_CACHE = {}
        return _REFERENCE_REACTIONS_CACHE
    
    try:
        with open(taxonomy_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        _REFERENCE_REACTIONS_CACHE = {}
        for entry in data:
            rxn_id = entry.get("id")
            refs = entry.get("reference_reactions", [])
            if rxn_id and refs:
                _REFERENCE_REACTIONS_CACHE[rxn_id] = refs
        
        logger.debug(f"Loaded reference reactions for {len(_REFERENCE_REACTIONS_CACHE)} reaction types")
        return _REFERENCE_REACTIONS_CACHE
        
    except Exception as e:
        logger.warning(f"Failed to load reference reactions: {e}")
        _REFERENCE_REACTIONS_CACHE = {}
        return _REFERENCE_REACTIONS_CACHE


# ============================================================================
# DRFP Fingerprint Functions
# ============================================================================

def _compute_drfp_fingerprint(reaction_smiles: str):
    """
    Compute DRFP fingerprint for a reaction SMILES.
    
    Returns:
        numpy array or None if computation fails
    """
    try:
        from drfp import DrfpEncoder
        import numpy as np
        
        encoder = DrfpEncoder()
        fps = encoder.encode([reaction_smiles])
        return np.array(fps[0]) if fps else None
    except ImportError:
        logger.debug("DRFP not available")
        return None
    except Exception as e:
        logger.debug(f"DRFP computation failed: {e}")
        return None


def _get_reference_drfp_index() -> Dict[str, List[Tuple[str, Any]]]:
    """
    Get or build the DRFP index for reference reactions.
    
    Returns:
        Dict mapping reaction_type_id to list of (ref_smiles, drfp_vector)
    """
    global _REFERENCE_DRFP_CACHE
    
    if _REFERENCE_DRFP_CACHE is not None:
        return _REFERENCE_DRFP_CACHE
    
    try:
        import numpy as np
    except ImportError:
        _REFERENCE_DRFP_CACHE = {}
        return _REFERENCE_DRFP_CACHE
    
    refs = _load_reference_reactions()
    _REFERENCE_DRFP_CACHE = {}
    
    for rxn_type, ref_smiles_list in refs.items():
        drfp_list = []
        for ref_smiles in ref_smiles_list:
            drfp = _compute_drfp_fingerprint(ref_smiles)
            if drfp is not None:
                drfp_list.append((ref_smiles, drfp))
        if drfp_list:
            _REFERENCE_DRFP_CACHE[rxn_type] = drfp_list
    
    logger.debug(f"Built DRFP index for {len(_REFERENCE_DRFP_CACHE)} reaction types")
    return _REFERENCE_DRFP_CACHE


def clear_caches():
    """Clear all module caches. Useful for testing or reloading taxonomy."""
    global _TAXONOMY_SMARTS_CACHE, _REFERENCE_REACTIONS_CACHE, _REFERENCE_DRFP_CACHE
    _TAXONOMY_SMARTS_CACHE = None
    _REFERENCE_REACTIONS_CACHE = None
    _REFERENCE_DRFP_CACHE = None


# ============================================================================
# DRFP Similarity Detection
# ============================================================================

def detect_by_drfp_similarity(
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
        
    Example:
        >>> matches = detect_by_drfp_similarity("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
        >>> matches[0]
        ('suzuki_miyaura', 'Suzuki-Miyaura Coupling', 0.95)
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


# ============================================================================
# Taxonomy SMARTS Loading
# ============================================================================

def _load_taxonomy_smarts() -> Dict[str, Any]:
    """
    Load SMARTS patterns from reaction_types.json taxonomy file.
    
    DEPRECATED: The 'smarts' field has been removed from reaction_types.json.
    This function will return an empty dict. Use data-driven detection instead.
    
    Returns:
        Dict mapping reaction_type_id to {smarts, name, category, aliases}
        (empty if no smarts fields exist in taxonomy)
    """
    global _TAXONOMY_SMARTS_CACHE
    
    if _TAXONOMY_SMARTS_CACHE is not None:
        return _TAXONOMY_SMARTS_CACHE
    
    import json
    
    taxonomy_path = _get_taxonomy_path()
    
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
        
        if not _TAXONOMY_SMARTS_CACHE:
            logger.debug("No SMARTS patterns found in taxonomy (field removed)")
        else:
            logger.debug(f"Loaded {len(_TAXONOMY_SMARTS_CACHE)} reaction SMARTS from taxonomy")
        return _TAXONOMY_SMARTS_CACHE
        
    except Exception as e:
        logger.warning(f"Failed to load taxonomy SMARTS: {e}")
        _TAXONOMY_SMARTS_CACHE = {}
        return _TAXONOMY_SMARTS_CACHE


def get_taxonomy_entry(reaction_type_id: str) -> Optional[Dict[str, Any]]:
    """
    Get taxonomy entry for a specific reaction type.
    
    Args:
        reaction_type_id: The reaction type ID (e.g., 'suzuki_miyaura')
        
    Returns:
        Dict with smarts, name, category, aliases, or None if not found
    """
    taxonomy = _load_taxonomy_smarts()
    return taxonomy.get(reaction_type_id)


# ============================================================================
# Reaction Running Validation
# ============================================================================

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
        reactant_orderings = []
        if len(reactant_mols) <= 4:
            for perm in permutations(range(len(reactant_mols))):
                if len(perm) >= num_templates:
                    reactant_orderings.append([reactant_mols[i] for i in perm[:num_templates]])
                if len(reactant_orderings) >= 10:  # Limit permutations
                    break
        else:
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
                        Chem.SanitizeMol(pred_mol)
                        smi = Chem.MolToSmiles(pred_mol, canonical=True)
                        predicted_smiles.add(smi)
                    except Exception:
                        continue
                
                if not predicted_smiles:
                    continue
                
                # Calculate Jaccard similarity
                intersection = predicted_smiles & actual_product_smiles
                union = predicted_smiles | actual_product_smiles
                
                if intersection:
                    match_score = len(intersection) / len(union)
                    best_match_score = max(best_match_score, match_score)
                    
                    if match_score >= 0.8:
                        return True, match_score * 0.15
        
        if best_match_score >= 0.5:
            return True, best_match_score * 0.15
        
        return False, 0.0
        
    except Exception as e:
        logger.debug(f"Reaction validation failed: {e}")
        return False, 0.0


# ============================================================================
# Taxonomy SMARTS Detection (Main Function)
# ============================================================================

def detect_by_taxonomy_smarts(reaction_smiles: str) -> List[Tuple[str, str, float]]:
    """
    Detect reaction type by matching against taxonomy SMARTS patterns.
    
    DEPRECATED: This function is deprecated as the `smarts` field has been
    removed from reaction_types.json. Use data-driven detection via
    `chemtools.taxonomy.data_driven_detection.detect_by_reactants()` instead,
    or use DRFP similarity matching via `detect_by_drfp_similarity()`.
    
    Args:
        reaction_smiles: Full reaction SMILES (reactants>>products)
        
    Returns:
        Empty list (smarts field removed from taxonomy).
        
    See Also:
        - detect_by_drfp_similarity: For similarity-based detection
        - chemtools.taxonomy.data_driven_detection.detect_by_reactants: For data-driven detection
    """
    import warnings
    warnings.warn(
        "detect_by_taxonomy_smarts is deprecated. The 'smarts' field has been removed "
        "from reaction_types.json. Use detect_by_drfp_similarity() or "
        "data_driven_detection.detect_by_reactants() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    
    # Load taxonomy - will return empty if no smarts fields exist
    taxonomy = _load_taxonomy_smarts()
    if not taxonomy:
        logger.debug("No SMARTS patterns found in taxonomy - returning empty")
        return []
    
    # Original implementation follows (for backwards compatibility if smarts exist)
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.warning("RDKit not available for taxonomy SMARTS matching")
        return []
    
    matches: List[Tuple[str, str, float, bool]] = []
    
    # Parse reaction SMILES
    parts = reaction_smiles.split(">>")
    if len(parts) != 2:
        return []
    
    reactants_smiles, products_smiles = parts
    
    # Parse reactants and products as molecules
    try:
        reactant_mols = []
        for smi in reactants_smiles.split("."):
            mol = Chem.MolFromSmiles(smi)
            if mol:
                reactant_mols.append(mol)
        
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
            rxn_pattern = AllChem.ReactionFromSmarts(smarts)
            if rxn_pattern is None:
                continue
            
            num_reactant_templates = rxn_pattern.GetNumReactantTemplates()
            num_product_templates = rxn_pattern.GetNumProductTemplates()
            
            if num_reactant_templates == 0:
                continue
            
            # === STEP 1: Match reactant templates ===
            template_matches = []
            for i in range(num_reactant_templates):
                template = rxn_pattern.GetReactantTemplate(i)
                matching_reactants = set()
                for j, mol in enumerate(reactant_mols):
                    if mol.HasSubstructMatch(template):
                        matching_reactants.add(j)
                template_matches.append(matching_reactants)
            
            if any(len(m) == 0 for m in template_matches):
                continue
            
            # Check valid assignment
            reactants_matched = False
            if num_reactant_templates == 1:
                reactants_matched = len(template_matches[0]) > 0
            elif num_reactant_templates == 2:
                for r1 in template_matches[0]:
                    for r2 in template_matches[1]:
                        if r1 != r2:
                            reactants_matched = True
                            break
                    if reactants_matched:
                        break
                if not reactants_matched and len(reactant_mols) == 2:
                    if 0 in template_matches[0] and 1 in template_matches[1]:
                        reactants_matched = True
                    elif 1 in template_matches[0] and 0 in template_matches[1]:
                        reactants_matched = True
            else:
                reactants_matched = all(len(m) > 0 for m in template_matches)
            
            if not reactants_matched:
                continue
            
            # === STEP 2: Match product templates ===
            products_matched = True
            
            if num_product_templates > 0 and product_mols:
                product_template_matches = []
                for i in range(num_product_templates):
                    template = rxn_pattern.GetProductTemplate(i)
                    matching_products = set()
                    for j, mol in enumerate(product_mols):
                        if mol.HasSubstructMatch(template):
                            matching_products.add(j)
                    product_template_matches.append(matching_products)
                
                if any(len(m) == 0 for m in product_template_matches):
                    products_matched = False
                else:
                    if num_product_templates == 1:
                        products_matched = len(product_template_matches[0]) > 0
                    elif num_product_templates == 2:
                        products_matched = False
                        for p1 in product_template_matches[0]:
                            for p2 in product_template_matches[1]:
                                if p1 != p2:
                                    products_matched = True
                                    break
                            if products_matched:
                                break
                        if not products_matched and len(product_mols) == 2:
                            if 0 in product_template_matches[0] and 1 in product_template_matches[1]:
                                products_matched = True
                            elif 1 in product_template_matches[0] and 0 in product_template_matches[1]:
                                products_matched = True
                    else:
                        products_matched = all(len(m) > 0 for m in product_template_matches)
            
            # === STEP 3: Calculate confidence ===
            if reactants_matched and products_matched:
                specificity = min(len(smarts) / 150.0, 1.0)
                
                if len(reactant_mols) == num_reactant_templates:
                    specificity = min(specificity + 0.1, 1.0)
                
                if num_product_templates > 0 and products_matched:
                    specificity = min(specificity + 0.05, 1.0)
                
                # === STEP 4: Validate by running reaction ===
                reaction_validated, validation_boost = _validate_reaction_by_running(
                    rxn_pattern, reactant_mols, product_mols
                )
                
                if reaction_validated:
                    specificity = min(specificity + validation_boost, 1.0)
                
                confidence = 0.75 + (0.20 * specificity)
                
                matches.append((rxn_id, info["name"], confidence, reaction_validated))
                
        except Exception as e:
            logger.debug(f"Failed to match SMARTS for {rxn_id}: {e}")
            continue
    
    # Sort by validation first, then confidence
    matches.sort(key=lambda x: (x[3], x[2]), reverse=True)
    
    # Return without validation flag
    return [(m[0], m[1], m[2]) for m in matches]
