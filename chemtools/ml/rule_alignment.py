"""
Rule-Alignment Scoring for ML Recommendations

This module implements a scoring system that reranks ML-based recommendations
based on their similarity to rule-based matching results. The alignment score
helps identify ML recommendations that are consistent with established chemistry
rules and precedents.

Key Features:
- Multi-component alignment scoring (solvent, metal, ligand, base, etc.)
- Configurable component weights
- Best-match finding between ML and rule recommendations
- Detailed alignment explanations
"""

from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, field
import warnings


@dataclass
class AlignmentScore:
    """Alignment score between an ML recommendation and a rule entry."""
    ml_rank: int
    rule_rank: Optional[int]
    total_score: float
    component_scores: Dict[str, float] = field(default_factory=dict)
    weights: Dict[str, float] = field(default_factory=dict)
    matched_rule_index: Optional[int] = None
    reasoning: List[str] = field(default_factory=list)


# Default weights for rule-alignment scoring
DEFAULT_ALIGNMENT_WEIGHTS = {
    'solvent': 0.20,
    'metal': 0.25,
    'ligand': 0.20,
    'base': 0.20,
    'additive': 0.05,
    'temperature': 0.05,
    'concentration': 0.03,
    'time': 0.02,
}


def normalize_reagent_identifier(identifier: str) -> str:
    """
    Normalize a reagent identifier (name or CAS) for comparison.
    
    Args:
        identifier: Reagent name, CAS number, or abbreviation
    
    Returns:
        Normalized identifier (lowercase, stripped)
    """
    if not identifier:
        return ""
    
    # Convert to string and lowercase
    norm = str(identifier).lower().strip()
    
    # Remove common prefixes/suffixes for better matching
    norm = norm.replace("(s)", "").replace("(l)", "").replace("(g)", "")
    norm = norm.replace("anhydrous", "").strip()
    
    return norm


def extract_metal_from_core(core: str) -> Optional[str]:
    """
    Extract metal symbol from a catalytic core string.
    
    Args:
        core: Catalytic core (e.g., "Pd/XPhos", "Cu/L1")
    
    Returns:
        Metal symbol (e.g., "Pd", "Cu") or None
    """
    if not core:
        return None
    
    # Common metals in catalysis
    metals = ['Pd', 'Cu', 'Ni', 'Fe', 'Co', 'Ru', 'Rh', 'Ir', 'Pt', 'Au', 'Ag']
    
    core_upper = core.upper()
    for metal in metals:
        if metal in core_upper:
            return metal
    
    return None


def extract_ligand_from_core(core: str) -> Optional[str]:
    """
    Extract ligand from a catalytic core string.
    
    Args:
        core: Catalytic core (e.g., "Pd/XPhos", "Cu/L1")
    
    Returns:
        Ligand name or None
    """
    if not core or '/' not in core:
        return None
    
    parts = core.split('/')
    if len(parts) >= 2:
        return parts[1].strip()
    
    return None


def find_chemical_in_list(
    target: str,
    chemical_list: List[Dict[str, Any]],
    role_filter: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Find a chemical in a list by matching name, CAS, or SMILES.
    
    Args:
        target: Target identifier (name, CAS, or SMILES)
        chemical_list: List of chemical dictionaries
        role_filter: Optional role to filter by
    
    Returns:
        Matching chemical dict or None
    """
    if not target or not chemical_list:
        return None
    
    target_norm = normalize_reagent_identifier(target)
    
    for chem in chemical_list:
        # Filter by role if specified
        if role_filter and chem.get('role') != role_filter:
            continue
        
        # Check various identifier fields
        for field in ['name', 'cas', 'smiles', 'abbreviation', 'uid']:
            value = chem.get(field)
            if value:
                value_norm = normalize_reagent_identifier(str(value))
                # Check exact match or substring match
                if value_norm == target_norm or target_norm in value_norm or value_norm in target_norm:
                    return chem
    
    return None


def compute_component_score(
    ml_value: Any,
    rule_value: Any,
    component_type: str
) -> float:
    """
    Compute alignment score for a single component.
    
    Args:
        ml_value: Value from ML recommendation
        rule_value: Value from rule entry
        component_type: Type of component ('solvent', 'metal', etc.)
    
    Returns:
        Component score in [0, 1]
    """
    # Handle missing values
    if ml_value is None and rule_value is None:
        return 0.5  # Neutral score if both missing
    if ml_value is None or rule_value is None:
        return 0.0  # No match if one is missing
    
    # Numeric components (temperature, time, concentration)
    if component_type in ['temperature', 'time', 'concentration']:
        try:
            ml_num = float(ml_value)
            rule_num = float(rule_value)
            
            # Avoid division by zero
            if rule_num == 0:
                return 1.0 if ml_num == 0 else 0.0
            
            # Compute relative difference
            rel_diff = abs(ml_num - rule_num) / abs(rule_num)
            
            # Convert to score: 1.0 for exact match, decreasing with difference
            # Allow 20% tolerance for full score
            if rel_diff <= 0.2:
                return 1.0
            elif rel_diff <= 0.5:
                return 1.0 - (rel_diff - 0.2) / 0.3
            else:
                return max(0.0, 1.0 - rel_diff)
        except (ValueError, TypeError):
            # Fall back to string comparison
            pass
    
    # String/categorical components
    ml_str = normalize_reagent_identifier(str(ml_value))
    rule_str = normalize_reagent_identifier(str(rule_value))
    
    # Exact match
    if ml_str == rule_str:
        return 1.0
    
    # Substring match (partial credit)
    if ml_str in rule_str or rule_str in ml_str:
        return 0.7
    
    # No match
    return 0.0


def compute_rule_alignment_score(
    ml_recommendation: Dict[str, Any],
    rule_entry: Dict[str, Any],
    weights: Optional[Dict[str, float]] = None
) -> Tuple[float, Dict[str, float]]:
    """
    Compute alignment score between an ML recommendation and a rule entry.
    
    Args:
        ml_recommendation: ML recommendation dictionary
        rule_entry: Rule-based entry dictionary
        weights: Component weights (uses defaults if None)
    
    Returns:
        Tuple of (total_score, component_scores)
    """
    if weights is None:
        weights = DEFAULT_ALIGNMENT_WEIGHTS.copy()
    
    component_scores = {}
    
    # Extract ML recommendation components
    ml_chemicals = ml_recommendation.get('chemicals', [])
    ml_conditions = ml_recommendation.get('conditions', {})
    ml_summary = ml_recommendation.get('summary', {})
    
    # Extract rule entry components
    rule_chemicals = rule_entry.get('chemicals', [])
    rule_conditions = rule_entry.get('conditions', {})
    
    # 1. Solvent alignment
    ml_solvent = find_chemical_in_list(
        ml_summary.get('solvent', {}).get('cas') or ml_summary.get('solvent', {}).get('name'),
        ml_chemicals,
        role_filter='solvent'
    )
    rule_solvent = None
    # Search for solvent in rule chemicals
    for chem in rule_chemicals:
        if chem.get('role') == 'solvent':
            rule_solvent = chem
            break
    
    if ml_solvent and rule_solvent:
        ml_solv_id = ml_solvent.get('cas') or ml_solvent.get('name')
        rule_solv_id = rule_solvent.get('cas') or rule_solvent.get('name')
        component_scores['solvent'] = compute_component_score(ml_solv_id, rule_solv_id, 'solvent')
    elif not ml_solvent and not rule_solvent:
        component_scores['solvent'] = 0.5  # Both missing - neutral
    else:
        component_scores['solvent'] = 0.0
    
    # 2. Metal alignment
    ml_core = ml_summary.get('core', '')
    rule_core = rule_entry.get('core') or ''
    
    ml_metal = extract_metal_from_core(ml_core)
    rule_metal = extract_metal_from_core(rule_core)
    
    component_scores['metal'] = compute_component_score(ml_metal, rule_metal, 'metal')
    
    # 3. Ligand alignment
    ml_ligand = extract_ligand_from_core(ml_core)
    rule_ligand = extract_ligand_from_core(rule_core)
    
    component_scores['ligand'] = compute_component_score(ml_ligand, rule_ligand, 'ligand')
    
    # 4. Base alignment
    ml_base = find_chemical_in_list(
        ml_summary.get('base', {}).get('cas') or ml_summary.get('base', {}).get('name'),
        ml_chemicals,
        role_filter='base'
    )
    rule_base = None
    # Search for base in rule chemicals
    for chem in rule_chemicals:
        if chem.get('role') == 'base':
            rule_base = chem
            break
    
    if ml_base and rule_base:
        ml_base_id = ml_base.get('cas') or ml_base.get('name')
        rule_base_id = rule_base.get('cas') or rule_base.get('name')
        component_scores['base'] = compute_component_score(ml_base_id, rule_base_id, 'base')
    elif not ml_base and not rule_base:
        component_scores['base'] = 0.5  # Both missing - neutral
    else:
        component_scores['base'] = 0.0
    
    # 5. Additive alignment
    ml_additive = find_chemical_in_list(
        None,
        ml_chemicals,
        role_filter='additive'
    )
    rule_additive = find_chemical_in_list(
        None,
        rule_chemicals,
        role_filter='additive'
    )
    
    if ml_additive and rule_additive:
        ml_add_id = ml_additive.get('cas') or ml_additive.get('name')
        rule_add_id = rule_additive.get('cas') or rule_additive.get('name')
        component_scores['additive'] = compute_component_score(ml_add_id, rule_add_id, 'additive')
    else:
        component_scores['additive'] = 0.5  # Neutral if both missing
    
    # 6. Temperature alignment
    ml_temp = None
    rule_temp = None
    
    if ml_conditions.get('temperature'):
        ml_temp = ml_conditions['temperature'].get('value')
    if rule_conditions.get('temperature'):
        rule_temp = rule_conditions['temperature'].get('value')
    
    component_scores['temperature'] = compute_component_score(ml_temp, rule_temp, 'temperature')
    
    # 7. Concentration alignment (if available)
    ml_conc = ml_conditions.get('concentration', {}).get('value')
    rule_conc = rule_conditions.get('concentration', {}).get('value')
    
    component_scores['concentration'] = compute_component_score(ml_conc, rule_conc, 'concentration')
    
    # 8. Time alignment
    ml_time = None
    rule_time = None
    
    if ml_conditions.get('time'):
        ml_time = ml_conditions['time'].get('value')
    if rule_conditions.get('time'):
        rule_time = rule_conditions['time'].get('value')
    
    component_scores['time'] = compute_component_score(ml_time, rule_time, 'time')
    
    # Compute weighted total score
    total_score = 0.0
    total_weight = 0.0
    
    for component, score in component_scores.items():
        weight = weights.get(component, 0.0)
        total_score += weight * score
        total_weight += weight
    
    # Normalize by total weight
    if total_weight > 0:
        total_score = total_score / total_weight
    
    return total_score, component_scores


def find_best_rule_match(
    ml_recommendation: Dict[str, Any],
    rule_recommendations: List[Dict[str, Any]],
    weights: Optional[Dict[str, float]] = None
) -> Tuple[int, float, Dict[str, float]]:
    """
    Find the best-matching rule entry for an ML recommendation.
    
    Args:
        ml_recommendation: ML recommendation dictionary
        rule_recommendations: List of rule-based recommendations
        weights: Component weights
    
    Returns:
        Tuple of (best_rule_index, best_score, component_scores)
    """
    if not rule_recommendations:
        return -1, 0.0, {}
    
    best_idx = -1
    best_score = -1.0
    best_components = {}
    
    for idx, rule_entry in enumerate(rule_recommendations):
        score, components = compute_rule_alignment_score(
            ml_recommendation,
            rule_entry,
            weights
        )
        
        if score > best_score:
            best_score = score
            best_idx = idx
            best_components = components
    
    return best_idx, best_score, best_components


def rerank_ml_by_rule_alignment(
    ml_recommendations: List[Dict[str, Any]],
    rule_recommendations: List[Dict[str, Any]],
    weights: Optional[Dict[str, float]] = None,
    ml_weight: float = 0.6,
    alignment_weight: float = 0.4
) -> List[Dict[str, Any]]:
    """
    Rerank ML recommendations based on alignment with rule-based results.
    
    Args:
        ml_recommendations: List of ML recommendations
        rule_recommendations: List of rule-based recommendations
        weights: Component weights for alignment scoring
        ml_weight: Weight for original ML score (default: 0.6)
        alignment_weight: Weight for rule-alignment score (default: 0.4)
    
    Returns:
        Reranked list of ML recommendations with alignment scores
    """
    if not ml_recommendations:
        return []
    
    if not rule_recommendations:
        warnings.warn("No rule recommendations available for alignment. Returning original ML ranking.")
        return ml_recommendations
    
    # Compute alignment scores for each ML recommendation
    scored_recommendations = []
    
    for ml_idx, ml_rec in enumerate(ml_recommendations):
        # Find best matching rule entry
        rule_idx, alignment_score, components = find_best_rule_match(
            ml_rec,
            rule_recommendations,
            weights
        )
        
        # Get original ML score (if available)
        original_score = 1.0 - (ml_idx / len(ml_recommendations))  # Decreasing score by rank
        if 'summary' in ml_rec and 'confidence' in ml_rec['summary']:
            original_score = ml_rec['summary']['confidence']
        
        # Compute combined score
        combined_score = ml_weight * original_score + alignment_weight * alignment_score
        
        # Build reasoning
        reasoning = []
        if alignment_score > 0.7:
            reasoning.append(f"Strong alignment with rule #{rule_idx + 1} (score: {alignment_score:.2f})")
        elif alignment_score > 0.5:
            reasoning.append(f"Moderate alignment with rule #{rule_idx + 1} (score: {alignment_score:.2f})")
        else:
            reasoning.append(f"Weak alignment with rules (best score: {alignment_score:.2f})")
        
        # Highlight strong component matches
        for comp, score in components.items():
            if score >= 0.9:
                reasoning.append(f"Excellent {comp} match ({score:.2f})")
        
        # Create alignment metadata
        alignment_meta = {
            'rule_alignment': {
                'best_match_index': rule_idx,
                'alignment_score': alignment_score,
                'component_scores': components,
                'weights': weights or DEFAULT_ALIGNMENT_WEIGHTS,
                'reasoning': reasoning
            },
            'combined_score': combined_score,
            'original_ml_score': original_score,
        }
        
        # Add alignment metadata to recommendation
        ml_rec_copy = ml_rec.copy()
        ml_rec_copy['alignment_meta'] = alignment_meta
        
        scored_recommendations.append((combined_score, ml_rec_copy))
    
    # Sort by combined score (descending)
    scored_recommendations.sort(key=lambda x: x[0], reverse=True)
    
    # Extract reranked recommendations and update ranks
    reranked = []
    for new_rank, (score, rec) in enumerate(scored_recommendations, start=1):
        rec['rank'] = new_rank
        rec['alignment_meta']['original_rank'] = rec.get('original_rank', rec.get('rank', new_rank))
        reranked.append(rec)
    
    return reranked


def explain_alignment(
    ml_recommendation: Dict[str, Any],
    rule_recommendation: Dict[str, Any],
    weights: Optional[Dict[str, float]] = None
) -> Dict[str, Any]:
    """
    Generate detailed explanation of alignment between ML and rule recommendations.
    
    Args:
        ml_recommendation: ML recommendation
        rule_recommendation: Rule recommendation
        weights: Component weights
    
    Returns:
        Explanation dictionary
    """
    total_score, component_scores = compute_rule_alignment_score(
        ml_recommendation,
        rule_recommendation,
        weights
    )
    
    explanation = {
        'total_alignment_score': total_score,
        'component_scores': component_scores,
        'weights': weights or DEFAULT_ALIGNMENT_WEIGHTS,
        'breakdown': []
    }
    
    # Generate detailed breakdown
    for component, score in component_scores.items():
        weight = (weights or DEFAULT_ALIGNMENT_WEIGHTS).get(component, 0.0)
        contribution = weight * score
        
        if score >= 0.9:
            match_quality = "Excellent"
        elif score >= 0.7:
            match_quality = "Good"
        elif score >= 0.5:
            match_quality = "Moderate"
        elif score >= 0.3:
            match_quality = "Weak"
        else:
            match_quality = "Poor"
        
        explanation['breakdown'].append({
            'component': component,
            'score': score,
            'weight': weight,
            'contribution': contribution,
            'match_quality': match_quality
        })
    
    # Sort by contribution
    explanation['breakdown'].sort(key=lambda x: x['contribution'], reverse=True)
    
    return explanation
