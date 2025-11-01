"""Mechanism-based similarity scoring for cross-family recommendations.

This module provides reaction mechanism similarity calculations to improve
cross-family precedent quality by considering catalytic systems and bond formations.
"""

from typing import Dict, Tuple


def calculate_mechanism_similarity(family_a: str, family_b: str) -> float:
    """Calculate mechanism similarity between two reaction families.
    
    Args:
        family_a: First reaction family (canonical name)
        family_b: Second reaction family (canonical name)
        
    Returns:
        Similarity score between 0.0 (unrelated) and 1.0 (identical)
    """
    if not family_a or not family_b:
        return 0.2  # Low similarity for unknown families
    
    if family_a == family_b:
        return 1.0
    
    # Normalize family names
    fam_a = family_a.lower().strip()
    fam_b = family_b.lower().strip()
    
    # Define mechanism similarity matrix based on:
    # - Catalyst system (Pd, Cu, Ni, etc.)
    # - Bond formation type (C-C, C-N, etc.)
    # - Reaction mechanism (oxidative addition, cross-coupling, etc.)
    MECHANISM_SIMILARITY = {
        # High similarity: Same catalyst family, similar mechanism
        ('suzuki', 'heck'): 0.85,  # Both Pd-catalyzed C-C coupling
        ('suzuki', 'negishi'): 0.80,  # Both Pd-catalyzed cross-coupling
        ('suzuki', 'stille'): 0.80,  # Both Pd-catalyzed cross-coupling
        ('heck', 'negishi'): 0.75,  # Both Pd-catalyzed
        ('heck', 'stille'): 0.75,  # Both Pd-catalyzed
        
        # Medium-high similarity: Same bond formation, different catalyst
        ('c_n_coupling_pd', 'c_n_coupling_cu'): 0.70,  # Same C-N bond, different metal
        ('c_n_coupling_pd', 'c_n_coupling_ni'): 0.65,  # Same C-N bond, different metal
        ('c_n_coupling_cu', 'c_n_coupling_ni'): 0.60,  # Same C-N bond, different metal
        
        # Medium similarity: Same catalyst, different bond formation
        ('c_n_coupling_pd', 'suzuki'): 0.65,  # Same Pd catalyst system
        ('c_n_coupling_pd', 'heck'): 0.60,  # Same Pd catalyst system
        ('c_n_coupling_pd', 'negishi'): 0.60,  # Same Pd catalyst system
        
        # Lower similarity: Different catalyst and mechanism but related
        ('suzuki', 'sonogashira'): 0.50,  # Both cross-coupling but different bonds
        ('c_n_coupling_cu', 'chan_lam'): 0.50,  # Both Cu-catalyzed C-N formation
        ('c_n_coupling_cu', 'ullmann'): 0.85,  # Ullmann is Cu-catalyzed C-N (high similarity)
        
        # Low similarity: Different mechanisms
        ('suzuki', 'amide_formation'): 0.30,  # Different mechanisms entirely
        ('heck', 'reductive_amination'): 0.25,  # Very different mechanisms
        ('c_n_coupling_pd', 'click_chemistry'): 0.20,  # Unrelated mechanisms
    }
    
    # Try both orderings of the family pair
    pair_ab = (fam_a, fam_b)
    pair_ba = (fam_b, fam_a)
    
    similarity = MECHANISM_SIMILARITY.get(pair_ab) or MECHANISM_SIMILARITY.get(pair_ba)
    
    if similarity is not None:
        return similarity
    
    # Fallback: Check for partial name matches (e.g., "suzuki" in "suzuki_miyaura")
    if _partial_family_match(fam_a, fam_b):
        return 0.40  # Moderate similarity for related families
    
    return 0.25  # Default low similarity for unrelated families


def _partial_family_match(family_a: str, family_b: str) -> bool:
    """Check if two family names have significant overlap."""
    # Split on common separators
    tokens_a = set(family_a.replace('_', ' ').replace('-', ' ').split())
    tokens_b = set(family_b.replace('_', ' ').replace('-', ' ').split())
    
    # Remove common words that don't indicate mechanism similarity
    common_words = {'coupling', 'formation', 'reaction', 'catalyzed'}
    tokens_a -= common_words
    tokens_b -= common_words
    
    if not tokens_a or not tokens_b:
        return False
    
    # Check for significant overlap (≥50% of tokens in common)
    overlap = len(tokens_a & tokens_b)
    min_tokens = min(len(tokens_a), len(tokens_b))
    
    return overlap / min_tokens >= 0.5


def enhance_precedent_similarity(
    precedent: Dict,
    detected_family: str,
    base_similarity: float,
    mechanism_weight: float = 0.3
) -> float:
    """Enhance precedent similarity score with mechanism awareness.
    
    Args:
        precedent: Precedent dictionary with rxn_type field
        detected_family: Query reaction's detected family
        base_similarity: Original DRFP/feature similarity score
        mechanism_weight: Weight for mechanism similarity (0.0-1.0)
        
    Returns:
        Enhanced similarity score combining base + mechanism similarity
    """
    if mechanism_weight <= 0:
        return base_similarity
    
    precedent_family = precedent.get('rxn_type', '')
    mechanism_sim = calculate_mechanism_similarity(detected_family, precedent_family)
    
    # Weighted combination: emphasize base similarity but boost with mechanism match
    enhanced_similarity = (
        (1 - mechanism_weight) * base_similarity + 
        mechanism_weight * mechanism_sim
    )
    
    # Ensure we don't exceed 1.0
    return min(enhanced_similarity, 1.0)


def get_compatibility_status(mechanism_similarity: float) -> str:
    """Get compatibility status based on mechanism similarity score.
    
    Args:
        mechanism_similarity: Similarity score between 0.0 and 1.0
        
    Returns:
        Compatibility status: 'compatible', 'moderate', or 'incompatible'
    """
    if mechanism_similarity >= 0.6:
        return 'compatible'
    elif mechanism_similarity >= 0.35:
        return 'moderate'
    else:
        return 'incompatible'


def get_representation_status(type_fraction: float, threshold: float = 0.15) -> str:
    """Get representation status based on reaction type fraction.
    
    Args:
        type_fraction: Fraction of precedents with this reaction type (0.0-1.0)
        threshold: Minimum threshold for well-represented status
        
    Returns:
        Representation status: 'well_represented' or 'underrepresented'
    """
    return 'well_represented' if type_fraction >= threshold else 'underrepresented'


def get_family_compatibility_score(families: list, primary_family: str) -> Dict[str, float]:
    """Get compatibility scores for a list of families relative to a primary family.
    
    Args:
        families: List of reaction families to score
        primary_family: Reference family for compatibility calculation
        
    Returns:
        Dictionary mapping family names to compatibility scores (0.0-1.0)
    """
    scores = {}
    for family in families:
        scores[family] = calculate_mechanism_similarity(primary_family, family)
    
    return scores