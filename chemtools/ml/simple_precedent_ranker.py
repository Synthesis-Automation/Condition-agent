"""
Simplified Precedent-Centric Recommendation

Strategy: Use precedents as primary source, then optionally rerank by
dataset analytics (popular reagents).

This is much simpler than the full fusion approach.
"""

from typing import Dict, List, Any, Optional, Tuple
from collections import Counter
import warnings


def rerank_by_analytics(
    precedents: List[dict],
    similarities: List[float],
    family: str
) -> Tuple[List[dict], List[float], List[str]]:
    """
    Rerank precedents by dataset analytics (popular reagents).
    
    Strategy:
    - Get most common catalysts, bases, solvents from dataset
    - Boost precedents using popular reagents
    
    Args:
        precedents: List of precedent reactions
        similarities: DRFP similarity scores (0-1)
        family: Reaction family name
    
    Returns:
        (reranked_precedents, reranked_scores, reasoning_messages)
    """
    from chemtools import chem
    
    reasoning = []
    
    try:
        # Get dataset statistics
        top_cores = chem.analytics.get_condition_cores(family, top_n=10)
        top_bases = chem.analytics.get_common_bases(family, top_n=10, min_yield=None)
        top_solvents = chem.analytics.get_common_solvents(family, top_n=10, min_yield=None)
        
        # Create popularity ranks (higher rank = more popular)
        core_ranks = {core: (len(top_cores) - i) for i, (core, _, _) in enumerate(top_cores)}
        base_ranks = {base: (len(top_bases) - i) for i, (base, _, _) in enumerate(top_bases)}
        solvent_ranks = {solv: (len(top_solvents) - i) for i, (solv, _, _) in enumerate(top_solvents)}
        
        max_rank = max(len(top_cores), len(top_bases), len(top_solvents))
        
        reasoning.append(f"Dataset analytics: {len(top_cores)} cores, {len(top_bases)} bases, {len(top_solvents)} solvents")
        
        # Score each precedent by reagent popularity
        boosted_scores = []
        for prec, sim in zip(precedents, similarities):
            boost = 0.0
            
            # Check catalyst popularity
            prec_core = prec.get('core', '')
            if prec_core in core_ranks:
                # Normalize to 0-0.3 range
                boost += (core_ranks[prec_core] / max_rank) * 0.3
            
            # Check base popularity
            prec_base = prec.get('base_uid', '')
            if prec_base in base_ranks:
                boost += (base_ranks[prec_base] / max_rank) * 0.2
            
            # Check solvent popularity
            prec_solvent = prec.get('solvent_uid', '')
            if prec_solvent in solvent_ranks:
                boost += (solvent_ranks[prec_solvent] / max_rank) * 0.2
            
            # Combine similarity with popularity boost
            combined_score = sim * (1.0 + boost)
            boosted_scores.append(combined_score)
        
        # Sort by combined score
        sorted_triples = sorted(
            zip(precedents, boosted_scores, similarities),
            key=lambda x: x[1],
            reverse=True
        )
        
        reranked_prec = [p for p, _, _ in sorted_triples]
        reranked_scores = [s for _, _, s in sorted_triples]  # Return original similarity scores
        
        reasoning.append(f"Reranked {len(precedents)} precedents by reagent popularity")
        
        return reranked_prec, reranked_scores, reasoning
        
    except Exception as e:
        warnings.warn(f"Analytics-based reranking failed: {e}")
        reasoning.append(f"Analytics reranking error: {e} - using similarity only")
        sorted_pairs = sorted(
            zip(precedents, similarities),
            key=lambda x: x[1],
            reverse=True
        )
        return (
            [p for p, _ in sorted_pairs],
            [s for _, s in sorted_pairs],
            reasoning
        )


def recommend_simple(
    reaction_smiles: str,
    family: str,
    k: int = 50,
    rerank_strategy: str = 'analytics',  # 'analytics' or 'none'
    relax: Optional[Dict] = None,
    filter_unknown_reagents: bool = False
) -> Dict[str, Any]:
    """
    Simple precedent-centric recommendation with optional reranking.
    
    Workflow:
    1. Find k similar precedents by DRFP similarity
    2. Optionally filter precedents with unknown reagents
    3. Optionally rerank by dataset analytics
    4. Extract top conditions from reranked precedents
    
    Args:
        reaction_smiles: Reaction SMILES string
        family: Reaction family name
        k: Number of precedents to retrieve
        rerank_strategy: 'analytics' (popular reagents) or 'none' (similarity only)
        relax: Precedent search parameters
        filter_unknown_reagents: If True, remove precedents with unknown base/solvent reagents
            Note: Only filters base and solvent, not catalyst cores (cores are complex "Metal/Ligand" strings)
    
    Returns:
        Recommendations with precedents and reasoning
    """
    from chemtools import precedent
    
    # Step 1: Get precedents by similarity
    relax = relax or {}
    relax.setdefault('reaction_smiles', reaction_smiles)
    relax.setdefault('use_drfp', True)
    relax.setdefault('precompute_drfp', False)
    
    pack = precedent.knn(family=family, features={}, k=k, relax=relax)
    precs = list(pack.get('precedents', []))
    
    # Extract similarity scores
    similarities = [
        p.get('drfp_similarity', p.get('similarity', 0.5)) 
        for p in precs
    ]
    
    reasoning = [f"Found {len(precs)} similar precedents (top similarity: {max(similarities) if similarities else 0:.3f})"]
    
    # Step 2: Filter unknown reagents (optional)
    if filter_unknown_reagents and precs:
        try:
            from chemtools import reagent
            
            filtered_precs = []
            filtered_sims = []
            removed_count = 0
            
            for prec, sim in zip(precs, similarities):
                # Check if all reagents in this precedent are in the database
                core = prec.get('core')
                base_uid = prec.get('base_uid')
                solvent_uid = prec.get('solvent_uid')
                
                # More lenient filtering: Only check base and solvent (not core)
                # Cores are often complex "Metal/Ligand" strings that won't be in the database
                all_known = True
                
                # Check base (if present)
                if base_uid:
                    result = reagent.enrich_reagent_info(str(base_uid), 'base')
                    if not result.get('found', False):
                        all_known = False
                
                # Check solvent (if present)
                if all_known and solvent_uid:
                    result = reagent.enrich_reagent_info(str(solvent_uid), 'solvent')
                    if not result.get('found', False):
                        all_known = False
                
                if all_known:
                    filtered_precs.append(prec)
                    filtered_sims.append(sim)
                else:
                    removed_count += 1
            
            precs = filtered_precs
            similarities = filtered_sims
            if removed_count > 0:
                reasoning.append(f"Filtered {removed_count} precedents with unknown base/solvent reagents (not in database)")
                
        except Exception as e:
            import warnings
            warnings.warn(f"Unknown reagent filtering failed: {e}. Using all precedents.")
            reasoning.append(f"Reagent filtering error: {e} - using all precedents")
    
    # Step 3: Rerank precedents (optional)
    if rerank_strategy == 'analytics':
        precs, similarities, rerank_reasons = rerank_by_analytics(
            precs, similarities, family
        )
        reasoning.extend(rerank_reasons)
    else:
        # No reranking, just sort by similarity
        sorted_pairs = sorted(
            zip(precs, similarities),
            key=lambda x: x[1],
            reverse=True
        )
        precs = [p for p, _ in sorted_pairs]
        similarities = [s for _, s in sorted_pairs]
        reasoning.append("Using similarity ranking only (no reranking)")
    
    # Step 4: Extract top conditions from precedents
    # Use Counter to find most common reagents
    cores = [p.get('core') for p in precs[:20] if p.get('core')]
    bases = [p.get('base_uid') for p in precs[:20] if p.get('base_uid')]
    solvents = [p.get('solvent_uid') for p in precs[:20] if p.get('solvent_uid')]
    
    top_cores = Counter(cores).most_common(5)
    top_bases = Counter(bases).most_common(5)
    top_solvents = Counter(solvents).most_common(5)
    
    reasoning.append(f"Top reagents from precedents: {len(top_cores)} cores, {len(top_bases)} bases, {len(top_solvents)} solvents")
    
    return {
        'precedents': precs,
        'similarities': similarities,
        'top_cores': top_cores,
        'top_bases': top_bases,
        'top_solvents': top_solvents,
        'reasoning': reasoning,
        'method': f'simple_precedent_ranking_{rerank_strategy}',
        'rerank_strategy': rerank_strategy
    }
