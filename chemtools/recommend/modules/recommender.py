"""
Main DRFP-based recommendation engine.

This module contains the primary recommendation function that uses reaction
fingerprinting (DRFP) for similarity-based precedent search.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional
from collections import Counter

from ...smiles import normalize_reaction
from ...router import detect_family
from ... import precedent, explain
from ..utils import canonical_family, median, pick_with_constraints

# Optional rxn-insight integration for ML-based family detection
try:
    from ...reaction_type_detector import detect_reaction_type as _rxn_detect
    _HAS_RXN_INSIGHT = True
except Exception:
    _HAS_RXN_INSIGHT = False


def recommend_from_reaction(
    reaction: str,
    k: int = 25,
    relax: Dict[str, Any] | None = None,
    constraint_rules: Dict[str, Any] | None = None,
    *,
    family_override: Optional[str] = None,
    max_variants: int = 3,
    rerank_strategy: str = 'rule',
    filter_unknown_reagents: bool = False,
    _build_formatted_output_fn=None,  # Injected to avoid circular imports
) -> Dict[str, Any]:
    """
    Recommend reaction conditions from a reaction SMILES string.
    
    This is the main recommendation function that uses DRFP-based reaction similarity.
    Precedents are ranked by similarity, with optional reranking by rules or analytics
    to prevent dataset quality issues.
    
    Args:
        reaction: Reaction SMILES string (reactants>>products or reactants>>products>>agents)
        k: Number of precedents to retrieve (default: 25)
        relax: Optional relaxation parameters for precedent search (passed to precedent.knn())
        constraint_rules: Optional constraint rules (inventory, blacklist, environmental, etc.)
        family_override: Override auto-detected reaction family
        max_variants: Maximum number of condition variants to generate (default: 3)
        rerank_strategy: Precedent reranking strategy (default: 'rule')
            - 'rule': Boost precedents matching chemical rules (e.g., Ullmann→Cu)
            - 'analytics': Boost precedents using popular reagents from dataset
            - 'none': Use DRFP similarity only (no reranking)
        filter_unknown_reagents: If True, remove precedents with unknown base/solvent reagents (default: False)
            Note: Only filters base and solvent, not catalyst cores (cores are complex "Metal/Ligand" strings)
    
    Returns:
        Dict with keys:
            - input: Normalized reaction SMILES
            - family: Detected reaction family
            - features: Feature dict (usually empty when DRFP is used)
            - recommendation: Primary recommendation (core, base, solvent, T_C, time_h, confidence)
            - alternatives: Alternative cores, bases, solvents
            - precedent_pack: Full precedent search results
            - reasons: Human-readable explanation
            - formatted: Formatted output with multiple variants
    
    Example:
        >>> from chemtools.recommend import recommend_from_reaction
        >>> results = recommend_from_reaction(
        ...     reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        ...     k=50,
        ...     rerank_strategy='rule',
        ...     filter_unknown_reagents=True
        ... )
        >>> print(results['recommendation']['core'])
        'Cu/phen'  # Correctly identifies Cu for Ullmann, not Pd
    """
    relax = dict(relax or {})
    max_variants = max(1, int(max_variants or 1))
    fam_override_clean = (family_override.strip() if family_override else None)

    # 1) Normalize reaction SMILES
    norm = normalize_reaction(reaction)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]

    # 2) Detect reaction family
    # Try rxn-insight ML-based detection first (if available), fallback to rule-based
    detection_source = "user_supplied" if fam_override_clean else "auto"
    fam = fam_override_clean or "Unknown"
    rxn_smiles_norm = norm.get("normalized") or reaction
    
    # Check if rxn-insight should be used
    use_rxn_insight = relax.get("use_rxn_insight")
    if use_rxn_insight is None:
        import os as _os
        env_off = (_os.environ.get("CHEMTOOLS_DISABLE_RXN_INSIGHT", "").strip().lower() 
                   in {"1", "true", "yes", "on"})
        use_rxn_insight = not env_off
    
    auto_family = None
    if bool(use_rxn_insight) and _HAS_RXN_INSIGHT:
        try:
            rxn_auto = _rxn_detect(rxn_smiles_norm)
            if rxn_auto and rxn_auto.get("success") and rxn_auto.get("mapped_family"):
                auto_family = str(rxn_auto.get("mapped_family") or "Unknown")
        except Exception:
            pass
    
    # Rule-based detection as fallback
    rule_info = detect_family(reactants)
    rule_family = rule_info.get("family") or "Unknown"
    
    # Prioritize: user override > ML detection > rule-based
    if not fam_override_clean:
        fam = auto_family or rule_family or "Unknown"
    else:
        fam = fam_override_clean
    
    fam = canonical_family(fam)
    
    # 3) Features: Keep empty for DRFP-based search (default)
    # DRFP uses reaction_smiles directly, no need for complex substrate featurization
    features: Dict[str, Any] = {}
    
    # 4) Retrieve precedents using DRFP-based k-NN search
    relax.setdefault("reaction_smiles", rxn_smiles_norm)
    relax.setdefault("use_drfp", True)  # Enable DRFP by default
    relax.setdefault("precompute_drfp", False)  # Use binary NPZ files (faster)
    relax.setdefault("selective_loading", True)  # Load only needed family
    
    pack = precedent.knn(family=fam, features=features, k=int(k), relax=relax)
    
    precs: List[Dict[str, Any]] = list(pack.get("precedents") or [])
    support = int(pack.get("support") or len(precs))
    
    # Initialize reasoning list for tracking filtering and reranking
    rerank_reasons = []
    
    # 5) Optional filtering: Remove precedents with unknown reagents
    if filter_unknown_reagents and precs:
        try:
            from ... import reagent
            
            filtered_precs = []
            removed_count = 0
            
            for prec in precs:
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
                else:
                    removed_count += 1
            
            precs = filtered_precs
            if removed_count > 0:
                rerank_reasons.append(f"Filtered {removed_count} precedents with unknown base/solvent reagents")
                
        except Exception as e:
            import warnings
            warnings.warn(f"Unknown reagent filtering failed: {e}. Using all precedents.")
            rerank_reasons.append(f"Reagent filtering error: {e} - using all precedents")
    
    # 6) Optional reranking: Boost precedents by rule match or reagent popularity
    if rerank_strategy in ['rule', 'analytics'] and precs:
        try:
            from ...ml.simple_precedent_ranker import rerank_by_rules, rerank_by_analytics
            
            # Extract similarity scores
            similarities = [
                p.get('drfp_similarity', p.get('similarity', 0.5)) 
                for p in precs
            ]
            
            if rerank_strategy == 'rule':
                precs, similarities, rerank_reasons = rerank_by_rules(
                    precs, similarities, rxn_smiles_norm, fam
                )
            elif rerank_strategy == 'analytics':
                precs, similarities, rerank_reasons = rerank_by_analytics(
                    precs, similarities, fam
                )
            
            # Update pack with reranked precedents
            pack['precedents'] = precs
            pack['similarities'] = similarities
            
        except Exception as e:
            import warnings
            warnings.warn(f"Precedent reranking failed: {e}. Using similarity only.")
            rerank_reasons = [f"Reranking error: {e} - using similarity only"]
    
    # 7) Vote for catalytic core (Laplace smoothing for better confidence)
    core_counts = Counter([str(p.get("core") or "") for p in precs if p.get("core")])
    labels = list(core_counts.keys())
    alpha = 1.0  # Laplace smoothing parameter
    denom = sum(core_counts.values()) + alpha * max(1, len(labels))
    scores = {L: (core_counts.get(L, 0) + alpha) / denom for L in labels}
    
    chosen_core = max(scores, key=scores.get) if scores else None
    core_vote_share = (
        (core_counts.get(chosen_core, 0) / max(1, sum(core_counts.values())))
        if chosen_core else 0.0
    )
    
    # 8) Choose base and solvent (conditioned on chosen core)
    if chosen_core:
        group = [p for p in precs if str(p.get("core") or "") == chosen_core]
    else:
        group = precs
    
    bases = [str(p.get("base_uid") or "") for p in group if p.get("base_uid")]
    solvents = [str(p.get("solvent_uid") or "") for p in group if p.get("solvent_uid")]
    
    base_counts = Counter(bases)
    solv_counts = Counter(solvents)
    
    base_list = [b for b, _ in base_counts.most_common()] or [
        str(p.get("base_uid") or "") for p in precs if p.get("base_uid")
    ]
    solv_list = [s for s, _ in solv_counts.most_common()] or [
        str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")
    ]
    
    # Apply constraints (inventory, blacklist, environmental, etc.)
    base_pick, base_filter = pick_with_constraints(base_list, constraint_rules or {})
    solv_pick, solv_filter = pick_with_constraints(solv_list, constraint_rules or {})
    
    # Fallback to all precedents if filtered out or empty
    if (not base_pick) and precs:
        all_bases = [str(p.get("base_uid") or "") for p in precs if p.get("base_uid")]
        base_pick, base_filter = pick_with_constraints(
            list(dict.fromkeys(all_bases)), constraint_rules or {}
        )
    if (not solv_pick) and precs:
        all_solv = [str(p.get("solvent_uid") or "") for p in precs if p.get("solvent_uid")]
        solv_pick, solv_filter = pick_with_constraints(
            list(dict.fromkeys(all_solv)), constraint_rules or {}
        )
    
    # 9) Compute median temperature and time from same-core group
    def nums(key: str, items: List[Dict[str, Any]]):
        return [p.get(key) for p in items if isinstance(p.get(key), (int, float))]
    
    T_med = median(nums("T_C", group) or nums("T_C", precs))
    t_med = median(nums("time_h", group) or nums("time_h", precs))
    
    # 10) Confidence score (based on core vote share)
    conf = 0.95 * core_vote_share if support >= 5 else 0.5 * core_vote_share
    conf = max(0.3, min(0.95, conf))
    
    # 11) Generate explanation
    reasons_pack = explain.for_pack(pack, features)
    
    # 12) Build recommendation dict
    recommendation = {
        "core": chosen_core,
        "base_uid": base_pick,
        "solvent_uid": solv_pick,
        "T_C": T_med,
        "time_h": t_med,
        "confidence": round(float(conf), 3),
    }
    
    alternatives = {
        "cores": core_counts.most_common(3),
        "bases": base_counts.most_common(3),
        "solvents": solv_counts.most_common(3),
    }
    
    detection_meta = {
        "auto_family": auto_family,
        "rule_family": rule_family,
        "source": detection_source,
    }
    
    # 13) Build formatted output with variants (use injected function to avoid circular imports)
    if _build_formatted_output_fn is None:
        # Import here to avoid circular dependency
        from .output_builder import build_formatted_output as _build_formatted_output_fn
    
    formatted = _build_formatted_output_fn(
        norm=norm,
        reaction=reaction,
        fam=fam,
        recommendation=recommendation,
        alternatives=alternatives,
        precs=precs,
        group=group,
        chosen_core=chosen_core,
        base_pick=base_pick,
        solv_pick=solv_pick,
        T_med=T_med,
        t_med=t_med,
        conf=conf,
        base_counts=base_counts,
        solv_counts=solv_counts,
        max_variants=max_variants,
    )
    
    return {
        "input": rxn_smiles_norm,
        "family": fam,
        "features": features,
        "bin": pack.get("prototype_id"),
        "recommendation": recommendation,
        "alternatives": alternatives,
        "precedent_pack": pack,
        "reasons": reasons_pack,
        "formatted": formatted,
        "detection": detection_meta,
        "constraint_filters": {
            "base": base_filter,
            "solvent": solv_filter,
        },
    }
