"""
Prototype: Multi-Source Evidence Fusion for ML Recommendations

This prototype demonstrates how to integrate:
- Precedent search (DRFP k-NN)
- Dataset analytics (frequency and yield statistics)  
- Rule-based matching (SMARTS schemes)
- ML yield prediction (DRFP-based models)

Key innovation: Adaptive weighting prevents over-reliance on potentially
biased top-k precedents by balancing with dataset-level statistics.
"""

from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, field
from collections import Counter
import warnings


@dataclass
class Candidate:
    """A candidate reaction condition."""
    core: str
    base: str
    solvent: str
    T_C: float = 100.0
    time_h: float = 12.0
    source: str = 'unknown'  # 'precedent', 'analytics', 'rules', 'hybrid'


@dataclass
class ScoredCandidate:
    """A candidate with multi-source scoring."""
    candidate: Candidate
    total_score: float
    component_scores: Dict[str, float] = field(default_factory=dict)
    weights: Dict[str, float] = field(default_factory=dict)
    confidence: str = 'MEDIUM'
    supporting_precedents: int = 0
    dataset_frequency: int = 0
    dataset_yield: Optional[float] = None
    rule_match: str = 'none'
    predicted_yield: Optional[float] = None
    reasoning: List[str] = field(default_factory=list)


def compute_diversity_score(precedents: List[dict]) -> float:
    """
    Measure diversity of precedent conditions.
    
    High diversity (>0.7): Precedents cover many different conditions
    Low diversity (<0.3): Precedents very similar (potential batch effect)
    
    Args:
        precedents: List of precedent reactions
    
    Returns:
        Diversity score in [0, 1]
    """
    if len(precedents) < 2:
        return 1.0
    
    distances = []
    for i in range(len(precedents)):
        for j in range(i + 1, len(precedents)):
            p1, p2 = precedents[i], precedents[j]
            
            dist = 0.0
            # Different catalyst: major difference
            if p1.get('core') != p2.get('core'):
                dist += 1.0
            # Different base: moderate difference
            if p1.get('base_uid') != p2.get('base_uid'):
                dist += 0.5
            # Different solvent: moderate difference
            if p1.get('solvent_uid') != p2.get('solvent_uid'):
                dist += 0.5
            # Temperature difference
            T1 = p1.get('T_C', 80)
            T2 = p2.get('T_C', 80)
            if T1 and T2 and abs(T1 - T2) > 20:
                dist += 0.3
            
            distances.append(dist)
    
    if not distances:
        return 1.0
    
    # Normalize by max possible distance
    max_dist = 2.3  # 1.0 + 0.5 + 0.5 + 0.3
    avg_dist = sum(distances) / len(distances)
    diversity = min(1.0, avg_dist / max_dist)
    
    return diversity


def deduplicate_candidates(candidates: List[Candidate]) -> List[Candidate]:
    """
    Remove duplicate candidates based on (core, base, solvent) tuple.
    
    Keeps the first occurrence of each unique condition.
    
    Args:
        candidates: List of candidate conditions
    
    Returns:
        Deduplicated list
    """
    seen = set()
    unique = []
    
    for candidate in candidates:
        # Create tuple key (case-insensitive)
        key = (
            (candidate.core or '').strip().lower(),
            (candidate.base or '').strip().lower(),
            (candidate.solvent or '').strip().lower()
        )
        
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    
    return unique


def _normalize_reagent_name(name: str) -> str:
    """
    Normalize reagent name for matching.
    
    Handles:
    - Case differences
    - Whitespace variations
    - Common abbreviations
    
    Args:
        name: Reagent name or CAS number
    
    Returns:
        Normalized string
    """
    if not name:
        return ''
    
    normalized = str(name).strip().lower()
    
    # Remove common prefixes/suffixes
    normalized = normalized.replace('cas:', '').replace('cas ', '').strip()
    
    return normalized


def _match_catalytic_system(candidate_core: str, system_str: str) -> bool:
    """
    Match candidate core against catalytic system string.
    
    Handles formats like:
    - "Pd(OAc)2 + XPhos"
    - "Pd/XPhos"
    - "CuI"
    
    Args:
        candidate_core: Candidate catalyst core
        system_str: Catalytic system string from analytics
    
    Returns:
        True if match found
    """
    if not candidate_core or not system_str:
        return False
    
    cand_norm = _normalize_reagent_name(candidate_core)
    sys_norm = _normalize_reagent_name(system_str)
    
    # Exact match
    if cand_norm == sys_norm:
        return True
    
    # Partial match: candidate appears in system
    if cand_norm in sys_norm or sys_norm in cand_norm:
        return True
    
    # Split system string and check components
    # e.g., "Pd(OAc)2 + XPhos" -> ["pd(oac)2", "xphos"]
    components = []
    for sep in [' + ', '+', '/', ' / ']:
        if sep in sys_norm:
            components = [c.strip() for c in sys_norm.split(sep)]
            break
    
    if components:
        for comp in components:
            if cand_norm == comp or cand_norm in comp or comp in cand_norm:
                return True
    
    return False


def collect_evidence(
    reaction_smiles: str,
    family: str,
    k: int = 50,
    relax: Optional[Dict] = None,
) -> Dict[str, Any]:
    """
    Collect evidence from all sources.
    
    Args:
        reaction_smiles: Reaction SMILES string
        family: Reaction family name
        k: Number of precedents to retrieve
        relax: Relaxation parameters for precedent search
    
    Returns:
        Evidence dictionary with precedents, analytics, rules, and ML info
    """
    from chemtools import chem, precedent
    
    evidence = {
        'precedents': {},
        'analytics': {},
        'rules': {},
        'ml': {}
    }
    
    # 1. Precedent Search (DRFP-based k-NN)
    try:
        relax = relax or {}
        relax.setdefault('reaction_smiles', reaction_smiles)
        relax.setdefault('use_drfp', True)
        relax.setdefault('precompute_drfp', False)
        
        pack = precedent.knn(family=family, features={}, k=k, relax=relax)
        precs = list(pack.get('precedents', []))
        
        # Extract DRFP similarity scores if available
        similarities = [
            p.get('drfp_similarity', p.get('similarity', 0.5)) 
            for p in precs
        ]
        
        diversity = compute_diversity_score(precs)
        
        evidence['precedents'] = {
            'reactions': precs,
            'similarities': similarities,
            'diversity_score': diversity,
            'coverage': len(precs),
            'avg_similarity': sum(similarities) / len(similarities) if similarities else 0.0
        }
    except Exception as e:
        warnings.warn(f"Precedent search failed: {e}")
        evidence['precedents'] = {
            'reactions': [],
            'similarities': [],
            'diversity_score': 0.0,
            'coverage': 0,
            'avg_similarity': 0.0
        }
    
    # 2. Dataset Analytics
    try:
        stats = chem.analytics.get_stats(family)
        
        evidence['analytics'] = {
            'catalytic_systems': chem.analytics.get_common_catalytic_systems(
                family, top_n=20, min_yield=75.0
            ),
            'bases': chem.analytics.get_common_bases(
                family, top_n=15, min_yield=None
            ),
            'solvents': chem.analytics.get_common_solvents(
                family, top_n=15, min_yield=None
            ),
            'condition_cores': chem.analytics.get_condition_cores(
                family, top_n=15
            ),
            'dataset_size': stats.get('total_reactions', 0),
            'yield_stats': stats.get('yield_stats', {})
        }
    except Exception as e:
        warnings.warn(f"Analytics collection failed: {e}")
        evidence['analytics'] = {
            'catalytic_systems': [],
            'bases': [],
            'solvents': [],
            'condition_cores': [],
            'dataset_size': 0,
            'yield_stats': {}
        }
    
    # 3. Rule-Based Matching (placeholder for now)
    evidence['rules'] = {
        'matched_schemes': [],
        'confidence': 0.0,
        'recommended': None
    }
    
    # 4. ML Model Info
    evidence['ml'] = {
        'model_available': False,
        'model_confidence': 0.0,
        'predictor': None
    }
    
    return evidence


def compute_adaptive_weights(evidence: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compute adaptive weights α, β, γ, δ based on evidence quality.
    
    Principles:
        - Strong precedents (k>50, diversity>0.5) → High α
        - Large dataset (>10,000 reactions) → High β
        - Exact rule match → High γ
        - High ML confidence → High δ
    
    Args:
        evidence: Evidence dictionary from collect_evidence()
    
    Returns:
        Dict with 'weights' and 'reasoning'
    """
    # Default weights
    weights = {
        'α': 0.35,  # Precedent weight
        'β': 0.30,  # Analytics weight
        'γ': 0.20,  # Rule weight
        'δ': 0.15   # ML weight
    }
    
    reasoning = []
    
    # Extract evidence metrics
    n_prec = evidence['precedents']['coverage']
    diversity = evidence['precedents']['diversity_score']
    avg_sim = evidence['precedents']['avg_similarity']
    dataset_size = evidence['analytics']['dataset_size']
    rule_conf = evidence['rules']['confidence']
    ml_available = evidence['ml']['model_available']
    
    # Adjust based on precedent quality
    if n_prec < 10:
        weights['α'] *= 0.5
        weights['β'] *= 1.5
        reasoning.append(f"Few precedents (n={n_prec}) → rely more on analytics")
    elif n_prec > 50:
        weights['α'] *= 1.2
        reasoning.append(f"Many precedents (n={n_prec}) → trust precedent patterns")
    
    if diversity < 0.3:
        weights['α'] *= 0.7
        reasoning.append(f"Low diversity ({diversity:.2f}) → precedents may be biased")
    elif diversity > 0.7:
        weights['α'] *= 1.1
        reasoning.append(f"High diversity ({diversity:.2f}) → precedents representative")
    
    if avg_sim < 0.6:
        weights['α'] *= 0.8
        reasoning.append(f"Low similarity ({avg_sim:.2f}) → precedents less relevant")
    
    # Adjust based on dataset size
    if dataset_size > 10000:
        weights['β'] *= 1.3
        reasoning.append(f"Large dataset (n={dataset_size:,}) → trust analytics")
    elif dataset_size < 1000:
        weights['β'] *= 0.7
        reasoning.append(f"Small dataset (n={dataset_size:,}) → analytics less reliable")
    
    # Adjust based on rule match
    if rule_conf > 0.9:
        weights['γ'] *= 2.0
        reasoning.append(f"Strong rule match (conf={rule_conf:.2f}) → trust chemistry")
    elif rule_conf > 0.7:
        weights['γ'] *= 1.3
        reasoning.append(f"Moderate rule match (conf={rule_conf:.2f})")
    else:
        weights['γ'] *= 0.5
        reasoning.append("No strong rule match → low rule weight")
    
    # Adjust based on ML availability
    if not ml_available:
        weights['δ'] = 0.0
        reasoning.append("No ML model available → δ=0")
    else:
        ml_conf = evidence['ml'].get('model_confidence', 0.5)
        if ml_conf > 0.85:
            weights['δ'] *= 1.2
            reasoning.append(f"High ML confidence ({ml_conf:.2f})")
        elif ml_conf < 0.6:
            weights['δ'] *= 0.6
            reasoning.append(f"Low ML confidence ({ml_conf:.2f})")
    
    # Normalize to sum to 1.0
    total = sum(weights.values())
    if total > 0:
        weights = {k: v / total for k, v in weights.items()}
    
    return {
        'weights': weights,
        'reasoning': reasoning
    }


def score_from_precedents(
    candidate: Candidate,
    precedent_evidence: Dict[str, Any]
) -> float:
    """
    Score candidate based on precedent support.
    
    Combines:
        - Core match frequency in precedents
        - Base match frequency in precedents
        - Solvent match frequency in precedents
        - Weighted by similarity scores
    
    Args:
        candidate: Candidate condition
        precedent_evidence: Precedent evidence dict
    
    Returns:
        Precedent score in [0, 1]
    """
    precs = precedent_evidence['reactions']
    sims = precedent_evidence['similarities']
    
    if not precs:
        return 0.0
    
    # Count matches weighted by similarity
    core_score = 0.0
    base_score = 0.0
    solv_score = 0.0
    total_weight = 0.0
    
    for prec, sim in zip(precs, sims):
        weight = sim  # Use similarity as weight
        total_weight += weight
        
        if prec.get('core') == candidate.core:
            core_score += weight
        if prec.get('base_uid') == candidate.base:
            base_score += weight
        if prec.get('solvent_uid') == candidate.solvent:
            solv_score += weight
    
    if total_weight == 0:
        return 0.0
    
    # Average of three components
    prec_score = (core_score + base_score + solv_score) / (3 * total_weight)
    
    return min(1.0, prec_score)


def score_from_analytics(
    candidate: Candidate,
    analytics_evidence: Dict[str, Any]
) -> float:
    """
    Score candidate based on dataset analytics.
    
    Combines:
        - Frequency percentile (how common is this condition?)
        - Yield percentile (how successful is this condition?)
    
    Uses improved matching that handles:
        - Catalytic system string parsing (e.g., "Pd(OAc)2 + XPhos")
        - CAS number vs name matching
        - Partial string matching
    
    Args:
        candidate: Candidate condition
        analytics_evidence: Analytics evidence dict
    
    Returns:
        Analytics score in [0, 1]
    """
    # Look up candidate in analytics data
    systems = analytics_evidence['catalytic_systems']
    bases = analytics_evidence['bases']
    solvents = analytics_evidence['solvents']
    dataset_size = analytics_evidence['dataset_size']
    
    if dataset_size == 0:
        return 0.0
    
    # Find matching catalytic system (with improved matching)
    system_freq = 0
    system_yield = None
    for sys, count, avg_yield in systems:
        if _match_catalytic_system(candidate.core, sys):
            system_freq = count
            system_yield = avg_yield
            break
    
    # Find matching base (with normalization)
    base_freq = 0
    base_yield = None
    cand_base_norm = _normalize_reagent_name(candidate.base)
    for base, count, avg_yield in bases:
        base_norm = _normalize_reagent_name(base)
        if cand_base_norm == base_norm or cand_base_norm in base_norm or base_norm in cand_base_norm:
            base_freq = count
            base_yield = avg_yield
            break
    
    # Find matching solvent (with normalization)
    solv_freq = 0
    solv_yield = None
    cand_solv_norm = _normalize_reagent_name(candidate.solvent)
    for solv, count, avg_yield in solvents:
        solv_norm = _normalize_reagent_name(solv)
        if cand_solv_norm == solv_norm or cand_solv_norm in solv_norm or solv_norm in cand_solv_norm:
            solv_freq = count
            solv_yield = avg_yield
            break
    
    # Compute frequency score (normalized by dataset size)
    max_freq = max(
        max([c for _, c, _ in systems], default=1),
        max([c for _, c, _ in bases], default=1),
        max([c for _, c, _ in solvents], default=1)
    )
    
    freq_score = (
        (system_freq / max_freq if system_freq else 0) +
        (base_freq / max_freq if base_freq else 0) +
        (solv_freq / max_freq if solv_freq else 0)
    ) / 3.0
    
    # Compute yield score (average of available yields)
    yields = [y for y in [system_yield, base_yield, solv_yield] if y is not None]
    yield_score = (sum(yields) / len(yields) / 100.0) if yields else 0.5
    
    # Combine: 60% frequency, 40% yield
    analytics_score = 0.6 * freq_score + 0.4 * yield_score
    
    return min(1.0, analytics_score)


def score_from_rules(
    candidate: Candidate,
    rule_evidence: Dict[str, Any]
) -> float:
    """
    Score candidate based on rule matching.
    
    Args:
        candidate: Candidate condition
        rule_evidence: Rule evidence dict
    
    Returns:
        Rule score in [0, 1]
    """
    # Placeholder - would integrate with rule-based system
    return rule_evidence.get('confidence', 0.0)


def score_from_ml(
    candidate: Candidate,
    ml_evidence: Dict[str, Any],
    reaction_smiles: str
) -> float:
    """
    Score candidate using ML yield prediction.
    
    Args:
        candidate: Candidate condition
        ml_evidence: ML evidence dict
        reaction_smiles: Reaction SMILES
    
    Returns:
        ML score in [0, 1] (predicted yield / 100)
    """
    if not ml_evidence.get('model_available'):
        return 0.5  # Neutral score if no model
    
    predictor = ml_evidence.get('predictor')
    if predictor is None:
        return 0.5
    
    try:
        import pandas as pd
        
        # Build test dataframe row
        test_row = {
            'reaction_smiles': reaction_smiles,
            'core': candidate.core,
            'base_uid': candidate.base,
            'solvent_uid': candidate.solvent,
            'T_C': candidate.T_C,
            'time_h': candidate.time_h,
        }
        
        test_df = pd.DataFrame([test_row])
        yield_pred = predictor.predict(test_df)[0]
        
        return min(1.0, max(0.0, float(yield_pred) / 100.0))
    except Exception as e:
        warnings.warn(f"ML prediction failed: {e}")
        return 0.5


def score_candidate(
    candidate: Candidate,
    evidence: Dict[str, Any],
    weights: Dict[str, float],
    reaction_smiles: str
) -> ScoredCandidate:
    """
    Score a candidate using multi-source evidence fusion.
    
    Args:
        candidate: Candidate condition
        evidence: Complete evidence dict
        weights: Adaptive weights (α, β, γ, δ)
        reaction_smiles: Reaction SMILES
    
    Returns:
        ScoredCandidate with total score and components
    """
    # Compute component scores
    ps = score_from_precedents(candidate, evidence['precedents'])
    as_ = score_from_analytics(candidate, evidence['analytics'])
    rs = score_from_rules(candidate, evidence['rules'])
    ms = score_from_ml(candidate, evidence['ml'], reaction_smiles)
    
    # Fuse scores
    total_score = (
        weights['α'] * ps +
        weights['β'] * as_ +
        weights['γ'] * rs +
        weights['δ'] * ms
    )
    
    # Determine confidence level
    if total_score > 0.7:
        confidence = 'HIGH'
    elif total_score > 0.5:
        confidence = 'MEDIUM'
    else:
        confidence = 'LOW'
    
    # Build reasoning
    reasoning = []
    
    if ps > 0.5:
        n_support = sum(
            1 for p in evidence['precedents']['reactions']
            if (p.get('core') == candidate.core or
                p.get('base_uid') == candidate.base or
                p.get('solvent_uid') == candidate.solvent)
        )
        reasoning.append(f"Precedent support: {n_support} similar reactions")
    
    if as_ > 0.5:
        reasoning.append("Strong dataset analytics support")
    
    if rs > 0.7:
        reasoning.append("Matches chemistry rules")
    
    if ms > 0.7:
        reasoning.append(f"ML predicts {ms*100:.1f}% yield")
    
    return ScoredCandidate(
        candidate=candidate,
        total_score=total_score,
        component_scores={'PS': ps, 'AS': as_, 'RS': rs, 'MS': ms},
        weights=weights,
        confidence=confidence,
        reasoning=reasoning
    )


def recommend_with_fusion(
    reaction_smiles: str,
    family: str,
    k: int = 50,
    top_n: int = 5,
    relax: Optional[Dict] = None
) -> Dict[str, Any]:
    """
    Main recommendation function with multi-source evidence fusion.
    
    Args:
        reaction_smiles: Reaction SMILES string
        family: Reaction family name
        k: Number of precedents to retrieve
        top_n: Number of top recommendations to return
        relax: Relaxation parameters
    
    Returns:
        Recommendations with scores, evidence, and explanations
    """
    # Step 1: Collect evidence from all sources
    evidence = collect_evidence(reaction_smiles, family, k, relax)
    
    # Step 2: Compute adaptive weights
    weight_info = compute_adaptive_weights(evidence)
    weights = weight_info['weights']
    
    # Step 3: Generate candidates (simplified for prototype)
    candidates = []
    
    # From precedents
    for prec in evidence['precedents']['reactions'][:10]:
        candidates.append(Candidate(
            core=prec.get('core', 'Unknown'),
            base=prec.get('base_uid', 'Unknown'),
            solvent=prec.get('solvent_uid', 'Unknown'),
            T_C=prec.get('T_C', 100.0),
            time_h=prec.get('time_h', 12.0),
            source='precedent'
        ))
    
    # From analytics (top catalytic systems + top bases/solvents)
    for system, _, _ in evidence['analytics']['catalytic_systems'][:3]:
        for base, _, _ in evidence['analytics']['bases'][:2]:
            for solvent, _, _ in evidence['analytics']['solvents'][:2]:
                candidates.append(Candidate(
                    core=system,
                    base=base,
                    solvent=solvent,
                    T_C=100.0,
                    time_h=12.0,
                    source='analytics'
                ))
    
    # Step 3.5: Deduplicate candidates
    candidates = deduplicate_candidates(candidates)
    
    # Step 4: Score all candidates
    scored = []
    for candidate in candidates:
        scored_cand = score_candidate(candidate, evidence, weights, reaction_smiles)
        scored.append(scored_cand)
    
    # Step 5: Rank and return top-N
    ranked = sorted(scored, key=lambda x: x.total_score, reverse=True)
    
    return {
        'recommended_conditions': ranked[:top_n],
        'evidence': evidence,  # Raw evidence for conversion/analysis
        'evidence_summary': {
            'precedents': f"{evidence['precedents']['coverage']} found, "
                         f"diversity={evidence['precedents']['diversity_score']:.2f}, "
                         f"avg_similarity={evidence['precedents']['avg_similarity']:.2f}",
            'analytics': f"Dataset: {evidence['analytics']['dataset_size']:,} reactions",
            'rules': f"Confidence: {evidence['rules']['confidence']:.2f}",
            'ml': 'Available' if evidence['ml']['model_available'] else 'Not available'
        },
        'adaptive_weights': weight_info,
        'reasoning': weight_info.get('reasoning', []),
        'method': 'multi_source_fusion'
    }
