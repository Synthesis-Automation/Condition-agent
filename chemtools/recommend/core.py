"""
Core recommendation engine using DRFP-based reaction similarity.

This module provides the main recommendation functions that rely primarily on
reaction fingerprinting (DRFP) for similarity-based precedent search. Complex
substrate featurization is optional and only used when DRFP is explicitly disabled.

Key Philosophy:
- Use reaction-level similarity (DRFP) as the primary method
- Avoid complex substrate analysis and family-specific featurization
- Keep the code general and maintainable
- Delegate to precedent.knn() for similarity search

Public Functions:
    - recommend_from_reaction(): Main recommendation function
    - recommend_conditions_structured(): Structured output with multiple variants
"""

from __future__ import annotations

from typing import Dict, Any, List, Tuple, Optional
from collections import Counter

from ..smiles import normalize_reaction
from ..router import detect_family
from .. import precedent, constraints, explain
from .utils import canonical_family, median, pick_with_constraints

# Optional rxn-insight integration for ML-based family detection
try:
    from ..reaction_type_detector import detect_reaction_type as _rxn_detect
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
            from .. import reagent_lookup
            
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
                    result = reagent_lookup.enrich_reagent_info(str(base_uid), 'base')
                    if not result.get('found', False):
                        all_known = False
                
                # Check solvent (if present)
                if all_known and solvent_uid:
                    result = reagent_lookup.enrich_reagent_info(str(solvent_uid), 'solvent')
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
            from ..ml.simple_precedent_ranker import rerank_by_rules, rerank_by_analytics
            
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
    
    # 6) Vote for catalytic core (Laplace smoothing for better confidence)
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
    
    # 6) Choose base and solvent (conditioned on chosen core)
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
    
    # 7) Compute median temperature and time from same-core group
    def nums(key: str, items: List[Dict[str, Any]]):
        return [p.get(key) for p in items if isinstance(p.get(key), (int, float))]
    
    T_med = median(nums("T_C", group) or nums("T_C", precs))
    t_med = median(nums("time_h", group) or nums("time_h", precs))
    
    # 8) Confidence score (based on core vote share)
    conf = 0.95 * core_vote_share if support >= 5 else 0.5 * core_vote_share
    conf = max(0.3, min(0.95, conf))
    
    # 9) Generate explanation
    reasons_pack = explain.for_pack(pack, features)
    
    # 10) Build recommendation dict
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
    
    # 11) Build formatted output with variants
    formatted = _build_formatted_output(
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


def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    *,
    k: int = 50,
    limit: int = 5,
    relax: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Recommend conditions with structured output format for API/UI consumers.
    
    This is a wrapper around recommend_from_reaction() that provides
    a more structured output format with proper metadata and precedent summaries.
    
    Args:
        reaction: Reaction SMILES string
        reaction_type: Optional reaction type override
        k: Number of precedents (default: 50 for better coverage)
        limit: Maximum number of recommendations to return (default: 5)
        relax: Relaxation parameters for precedent search
        constraints: Constraint rules (inventory, blacklist, etc.)
    
    Returns:
        Dict with keys:
            - meta: Metadata (status, timestamp, strategy, result_count)
            - input: Normalized reaction SMILES and family info
            - detection: Family detection info with auto/rule sources
            - recommendations: List of condition variants (limited to `limit`)
            - alternatives: Alternative cores, bases, solvents
            - precedents: Precedent summary with top examples
    """
    import time
    from datetime import datetime
    
    # Start timing
    start_time = time.time()
    
    limit = max(1, int(limit or 1))
    cfg_relax = dict(relax or {})
    
    result = recommend_from_reaction(
        reaction=reaction,
        k=int(k or 50),
        relax=cfg_relax,
        constraint_rules=constraints or {},
        family_override=reaction_type,
        max_variants=limit,
    )
    
    # Calculate processing time
    processing_time_ms = round((time.time() - start_time) * 1000, 2)
    
    # Extract formatted section
    formatted = dict(result.get("formatted") or {})
    recommendations = list(formatted.get("recommended_conditions") or [])
    recommendations = recommendations[:limit]
    
    # Ensure rank is set correctly
    for idx, rec in enumerate(recommendations, start=1):
        rec.setdefault("rank", idx)
        summary = rec.get("summary")
        if isinstance(summary, dict):
            summary.setdefault("rank", idx)
    
    # Build detection section with all fields
    detection = dict(formatted.get("detection") or {})
    detected_family = detection.get("family") or result.get("family") or "Unknown"
    detection_confidence = detection.get("confidence", 0.95)
    
    if reaction_type and not detection.get("source"):
        detection["source"] = "user_supplied"
    detection.setdefault("source", detection.get("source") or "auto")
    detection["detected_reaction_type"] = detected_family
    detection["confidence"] = detection_confidence
    detection.setdefault("method", "drfp_precedent_search")
    detection.setdefault("provided_reaction_type", reaction_type)
    
    # Build meta section with all fields
    meta = dict(formatted.get("meta") or {})
    meta["generated_at"] = datetime.utcnow().isoformat() + "Z"
    meta.setdefault("model", "drfp_similarity")
    meta.setdefault("status", "success")
    meta.setdefault("strategy", "drfp_similarity")
    meta["result_count"] = len(recommendations)
    meta["processing_time_ms"] = processing_time_ms
    
    # Build input section with all fields
    input_data = dict(formatted.get("input") or {})
    input_data["reaction_smiles"] = input_data.get("reaction_smiles") or reaction
    input_data["requested_reaction_type"] = reaction_type
    input_data.setdefault("family", detected_family)
    
    # Build precedent summary
    pack = result.get("precedent_pack") or {}
    precedents = list(pack.get("precedents") or [])
    top_precedents = [
        {
            "reaction_id": p.get("reaction_id"),
            "core": p.get("core"),
            "yield_pct": p.get("yield_pct"),
        }
        for p in precedents[:10]
        if p.get("reaction_id")
    ]
    
    core_family = result.get("family")
    core_support = len([p for p in precedents if p.get("core") == core_family])
    
    precedent_summary = {
        "total_considered": len(precedents),
        "core_family": core_family,
        "core_support": core_support,
        "top_precedents": top_precedents,
    }
    
    # Include detailed precedent information from formatted output
    precedents_used = formatted.get("precedents_used")
    
    return {
        "meta": meta,
        "input": input_data,
        "detection": detection,
        "recommendations": recommendations,
        "alternatives": result.get("alternatives"),
        "precedents": precedent_summary,
        "precedents_used": precedents_used,
        "constraint_filters": result.get("constraint_filters"),
    }




def _convert_fusion_to_core_format(
    fusion_results: Dict[str, Any],
    reaction: str,
    rxn_smiles_norm: str,
    fam: str,
    features: Dict[str, Any],
    detection_meta: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Convert fusion recommender output to core.py format for API compatibility.
    
    Fusion format has:
        - recommended_conditions: List of scored candidates (ScoredCandidate dataclass objects)
        - evidence: Multi-source evidence data
        - adaptive_weights: Dynamic weight adjustments
        - reasoning: Human-readable explanations
    
    Core format has:
        - recommendation: Primary recommendation
        - alternatives: Alternative cores/bases/solvents
        - formatted: Multi-variant output
        - precedent_pack: Precedent search results
        - reasons: Explanations
    """
    from collections import Counter
    
    recommended_conditions_raw = fusion_results.get('recommended_conditions', [])
    evidence = fusion_results.get('evidence', {})
    weights = fusion_results.get('adaptive_weights', {}).get('weights', {})
    reasoning = fusion_results.get('reasoning', [])
    
    # Convert ScoredCandidate dataclass objects to dicts
    recommended_conditions = []
    for scored_cand in recommended_conditions_raw:
        # Check if it's a dataclass object (has __dataclass_fields__)
        if hasattr(scored_cand, '__dataclass_fields__'):
            # It's a dataclass - extract fields
            cand = scored_cand.candidate
            recommended_conditions.append({
                'candidate': {
                    'core': cand.core,
                    'base': cand.base,
                    'solvent': cand.solvent,
                    'temperature': cand.T_C,
                    'time': cand.time_h,
                },
                'total_score': scored_cand.total_score,
                'component_scores': scored_cand.component_scores,
                'confidence_level': scored_cand.confidence,
            })
        else:
            # Already a dict - use as-is
            recommended_conditions.append(scored_cand)
    
    # Extract precedent pack from evidence
    prec_evidence = evidence.get('precedents', {})
    precedents = prec_evidence.get('reactions', [])
    
    precedent_pack = {
        'precedents': precedents,
        'support': prec_evidence.get('coverage', len(precedents)),
        'similarities': prec_evidence.get('similarities', []),
        'diversity_score': prec_evidence.get('diversity_score'),
        'avg_similarity': prec_evidence.get('avg_similarity'),
        'method': 'fusion',
    }
    
    # Build primary recommendation from top result
    if recommended_conditions:
        top = recommended_conditions[0]
        candidate = top.get('candidate', {})
        
        recommendation = {
            'core': candidate.get('core'),
            'base_uid': candidate.get('base'),
            'solvent_uid': candidate.get('solvent'),
            'T_C': candidate.get('temperature'),
            'time_h': candidate.get('time'),
            'confidence': top.get('confidence_level', 0.5),
            'fusion_score': top.get('total_score'),
            'fusion_method': 'multi_source_evidence',
        }
    else:
        recommendation = {
            'core': None,
            'base_uid': None,
            'solvent_uid': None,
            'T_C': None,
            'time_h': None,
            'confidence': 0.0,
            'fusion_score': 0.0,
            'fusion_method': 'multi_source_evidence',
        }
    
    # Build alternatives from all recommendations (NOW using converted dicts)
    core_counts = Counter([
        c.get('candidate', {}).get('core')
        for c in recommended_conditions
        if c.get('candidate', {}).get('core')
    ])
    base_counts = Counter([
        c.get('candidate', {}).get('base')
        for c in recommended_conditions
        if c.get('candidate', {}).get('base')
    ])
    solv_counts = Counter([
        c.get('candidate', {}).get('solvent')
        for c in recommended_conditions
        if c.get('candidate', {}).get('solvent')
    ])
    
    alternatives = {
        'cores': core_counts.most_common(3),
        'bases': base_counts.most_common(3),
        'solvents': solv_counts.most_common(3),
    }
    
    # Build formatted output compatible with existing structure
    # Pass the converted dict version, not the raw fusion_results
    fusion_results_converted = {
        'recommended_conditions': recommended_conditions,  # Use converted version
        'evidence': evidence,
        'adaptive_weights': fusion_results.get('adaptive_weights', {}),
        'reasoning': reasoning,
    }
    
    formatted = _build_formatted_output_from_fusion(
        fusion_results=fusion_results_converted,
        reaction=reaction,
        rxn_smiles_norm=rxn_smiles_norm,
        fam=fam,
        recommendation=recommendation,
        alternatives=alternatives,
        precedents=precedents,
        weights=weights,
    )
    
    # Build reasons/explanations
    reasons_text = '\n'.join(reasoning) if reasoning else "Multi-source evidence fusion"
    
    return {
        'input': rxn_smiles_norm,
        'family': fam,
        'features': features,
        'bin': None,
        'recommendation': recommendation,
        'alternatives': alternatives,
        'precedent_pack': precedent_pack,
        'reasons': reasons_text,
        'formatted': formatted,
        'detection': detection_meta,
        'fusion_meta': {
            'adaptive_weights': weights,
            'reasoning': reasoning,
            'evidence_summary': {
                'precedents': prec_evidence.get('coverage', 0),
                'diversity': prec_evidence.get('diversity_score'),
                'dataset_size': evidence.get('analytics', {}).get('dataset_size', 0),
            },
        },
        'constraint_filters': {},
    }


def _build_formatted_output_from_fusion(
    fusion_results: Dict[str, Any],
    reaction: str,
    rxn_smiles_norm: str,
    fam: str,
    recommendation: Dict[str, Any],
    alternatives: Dict[str, Any],
    precedents: List[Dict[str, Any]],
    weights: Dict[str, float],
) -> Dict[str, Any]:
    """
    Build formatted output structure from fusion results.
    
    Converts fusion recommendations into the standard multi-variant format
    expected by the API.
    """
    from ..smiles import normalize_reaction
    
    norm = normalize_reaction(reaction)
    recommended_conditions = fusion_results.get('recommended_conditions', [])
    
    # Helper: Lookup reagent info
    def _lookup(uid: str, role: str) -> Dict[str, Any]:
        """Look up reagent information by CAS number or name."""
        try:
            from .. import reagent_lookup
            
            role_to_type = {
                "base": "base",
                "solvent": "solvent",
                "ligand": "ligand",
                "metal_precursor": "metal_precursor",
                "catalyst": "metal_precursor",
            }
            reagent_type = role_to_type.get(role, role)
            
            res = reagent_lookup.enrich_reagent_info(uid, reagent_type)
            if res and res.get("found"):
                return {
                    "name": res.get("name"),
                    "token": res.get("abbreviation"),
                    "uid": res.get("cas") or uid,
                    "smiles": res.get("smiles"),
                    "inchi_key": res.get("inchi_key"),
                }
        except Exception:
            pass
        
        return {
            "uid": uid,
            "name": f"[{role.replace('_', ' ').title()}] {uid}",
            "unknown": True
        }
    
    def _chemical_payload(uid: Optional[str], role: str) -> Dict[str, Any] | None:
        if not uid:
            return None
        rec = _lookup(uid, role)
        
        display_name = rec.get("name")
        if rec.get("unknown"):
            display_name = f"[Unknown {role.replace('_', ' ')}] CAS {rec.get('uid')}"
        
        return {
            "name": display_name or rec.get("token") or rec.get("uid") or uid,
            "abbreviation": rec.get("token") or None,
            "cas": rec.get("uid") or uid,
            "smiles": rec.get("smiles"),
            "equivalents": None,
            "role": role,
        }
    
    # Build reactant chemicals
    reactants_chems: List[Dict[str, Any]] = []
    for r in (norm.get("reactants") or []):
        smi = r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or ""
        reactants_chems.append({
            "name": None,
            "cas": None,
            "smiles": smi or None,
            "equivalents": None,
            "role": "starting_material",
        })
    
    # Convert each fusion recommendation to formatted variant
    variants: List[Dict[str, Any]] = []
    
    for idx, rec in enumerate(recommended_conditions, start=1):
        candidate = rec.get('candidate', {})
        scores = rec.get('component_scores', {})
        
        # Build chemical list
        chems = list(reactants_chems)
        
        # Add core (catalyst system)
        core = candidate.get('core')
        if core:
            core_payload = _chemical_payload(core, 'catalyst')
            if core_payload:
                chems.append(core_payload)
        
        # Add base
        base = candidate.get('base')
        base_payload = _chemical_payload(base, 'base')
        if base_payload:
            chems.append(base_payload)
        
        # Add solvent
        solvent = candidate.get('solvent')
        solvent_payload = _chemical_payload(solvent, 'solvent')
        if solvent_payload:
            chems.append(solvent_payload)
        
        # Conditions
        temp = candidate.get('temperature')
        time = candidate.get('time')
        cond_text = {
            "temperature": (f"{int(temp)} °C" if temp is not None else "80 °C"),
            "time": (f"{time} h" if time is not None else "12 h"),
            "atmosphere": None,
        }
        
        # Build summary with fusion metadata
        summary = {
            "rank": idx,
            "core": core,
            "base": base_payload,
            "solvent": solvent_payload,
            "confidence": rec.get('confidence_level', 0.5),
            "fusion_score": rec.get('total_score', 0.0),
            "component_scores": {
                "precedent": scores.get('PS', 0.0),
                "analytics": scores.get('AS', 0.0),
                "rules": scores.get('RS', 0.0),
                "ml": scores.get('MS', 0.0),
            },
            "adaptive_weights": weights,
            "support": {
                "count": len([p for p in precedents if p.get('core') == core]),
                "fraction_core": 0.0,
                "reference_population": len(precedents),
            },
        }
        
        variant = {
            "rank": idx,
            "reaction": {"smiles": rxn_smiles_norm},
            "chemicals": chems,
            "conditions": cond_text,
            "summary": summary,
            "combo": {
                "base_uid": base,
                "solvent_uid": solvent
            },
        }
        variants.append(variant)
    
    return {
        "meta": {
            "status": "success",
            "model": "fusion",
            "version": "2.0",
            "method": "multi_source_evidence_fusion",
            "adaptive_weights": weights,
        },
        "input": {
            "reaction_smiles": rxn_smiles_norm,
            "family": fam,
        },
        "detection": {
            "family": fam,
            "confidence": recommendation.get('confidence', 0.5),
        },
        "recommended_conditions": variants,
        "precedent_summary": {
            "total_precedents": len(precedents),
            "chosen_core": recommendation.get('core'),
            "core_support": len([p for p in precedents if p.get('core') == recommendation.get('core')]),
        },
    }


def _build_precedent_details(
    precs: List[Dict[str, Any]],
    chosen_core: Optional[str],
    group: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Build detailed precedent information for the recommendation output.
    
    Args:
        precs: All precedents from the search
        chosen_core: The chosen catalyst core
        group: Precedents matching the chosen core
        
    Returns:
        Dictionary with comprehensive precedent details
    """
    # Get top 10 precedents overall
    top_precedents = []
    for i, p in enumerate(precs[:10], 1):
        precedent_info = {
            "rank": i,
            "reaction_id": p.get("reaction_id"),
            "reaction_smiles": p.get("reaction_smiles"),
            "core": p.get("core") or p.get("condition_core"),
            "yield": p.get("yield"),
        }
        
        # Add detailed reagent information
        catalytic_system = p.get("catalytic_system", [])
        if catalytic_system:
            precedent_info["catalysts"] = [
                {
                    "name": cat.get("name"),
                    "cas": cat.get("cas"),
                    "role": "catalyst/ligand"
                }
                for cat in catalytic_system
            ]
        
        # Add reagents (bases, additives, etc.)
        reagents = p.get("reagents", [])
        if reagents:
            precedent_info["reagents"] = [
                {
                    "name": r.get("name"),
                    "cas": r.get("cas"),
                    "role": r.get("role", "reagent").lower()
                }
                for r in reagents
            ]
        
        # Add solvents
        solvents = p.get("solvents", [])
        if solvents:
            precedent_info["solvents"] = [
                {
                    "name": s.get("name"),
                    "cas": s.get("cas"),
                }
                for s in solvents
            ]
        
        # Add conditions
        conditions_data = {}
        if p.get("T_C") is not None:
            conditions_data["temperature_C"] = p.get("T_C")
        if p.get("time_h") is not None:
            conditions_data["time_h"] = p.get("time_h")
        
        # Add from conditions dict if available
        cond_dict = p.get("conditions", {})
        if isinstance(cond_dict, dict):
            if "temperature_c" in cond_dict:
                conditions_data["temperature_C"] = cond_dict["temperature_c"]
            if "time_h" in cond_dict:
                conditions_data["time_h"] = cond_dict["time_h"]
            if "yield_pct" in cond_dict:
                conditions_data["yield_pct"] = cond_dict["yield_pct"]
        
        if conditions_data:
            precedent_info["conditions"] = conditions_data
        
        # Add reference if available
        if p.get("reference"):
            precedent_info["reference"] = p.get("reference")
        
        top_precedents.append(precedent_info)
    
    # Get precedents matching the chosen core
    core_precedents = []
    if chosen_core and group:
        for i, p in enumerate(group[:5], 1):
            core_prec = {
                "rank": i,
                "reaction_id": p.get("reaction_id"),
                "core": p.get("core") or p.get("condition_core"),
                "base": p.get("base_uid"),
                "solvent": p.get("solvent_uid"),
                "yield": p.get("yield"),
            }
            
            # Add temperature and time if available
            if p.get("T_C") is not None:
                core_prec["temperature_C"] = p.get("T_C")
            if p.get("time_h") is not None:
                core_prec["time_h"] = p.get("time_h")
            
            core_precedents.append(core_prec)
    
    return {
        "total_count": len(precs),
        "top_precedents": top_precedents,
        "core_matched_precedents": {
            "core": chosen_core,
            "count": len(group) if group else 0,
            "examples": core_precedents,
        },
        "statistics": {
            "average_yield": _calculate_average_yield(precs),
            "yield_range": _calculate_yield_range(precs),
            "temperature_range": _calculate_temp_range(precs),
            "time_range": _calculate_time_range(precs),
        }
    }


def _calculate_average_yield(precs: List[Dict[str, Any]]) -> Optional[float]:
    """Calculate average yield from precedents."""
    yields = [p.get("yield") for p in precs if p.get("yield") is not None and isinstance(p.get("yield"), (int, float))]
    if not yields:
        return None
    return round(sum(yields) / len(yields), 1)


def _calculate_yield_range(precs: List[Dict[str, Any]]) -> Optional[List[float]]:
    """Calculate yield range from precedents."""
    yields = [p.get("yield") for p in precs if p.get("yield") is not None and isinstance(p.get("yield"), (int, float))]
    if not yields:
        return None
    return [round(min(yields), 1), round(max(yields), 1)]


def _calculate_temp_range(precs: List[Dict[str, Any]]) -> Optional[List[float]]:
    """Calculate temperature range from precedents."""
    temps = [p.get("T_C") for p in precs if p.get("T_C") is not None and isinstance(p.get("T_C"), (int, float))]
    if not temps:
        return None
    return [round(min(temps), 1), round(max(temps), 1)]


def _calculate_time_range(precs: List[Dict[str, Any]]) -> Optional[List[float]]:
    """Calculate time range from precedents."""
    times = [p.get("time_h") for p in precs if p.get("time_h") is not None and isinstance(p.get("time_h"), (int, float))]
    if not times:
        return None
    return [round(min(times), 1), round(max(times), 1)]


def _build_formatted_output(
    norm: Dict[str, Any],
    reaction: str,
    fam: str,
    recommendation: Dict[str, Any],
    alternatives: Dict[str, Any],
    precs: List[Dict[str, Any]],
    group: List[Dict[str, Any]],
    chosen_core: Optional[str],
    base_pick: Optional[str],
    solv_pick: Optional[str],
    T_med: Optional[float],
    t_med: Optional[float],
    conf: float,
    base_counts: Counter,
    solv_counts: Counter,
    max_variants: int,
) -> Dict[str, Any]:
    """
    Build formatted output with multiple condition variants.
    
    Internal helper function for generating API-friendly output format.
    """
    # Helper: Lookup reagent info from reagent database
    def _lookup(uid: str, role: str) -> Dict[str, Any]:
        """Look up reagent information by CAS number or name."""
        try:
            from .. import reagent_lookup
            
            # Map role to reagent type for lookup
            role_to_type = {
                "base": "base",
                "solvent": "solvent",
                "ligand": "ligand",
                "metal_precursor": "metal_precursor",
                "catalyst": "metal_precursor",
            }
            reagent_type = role_to_type.get(role, role)
            
            # Try to enrich the reagent info
            res = reagent_lookup.enrich_reagent_info(uid, reagent_type)
            if res and res.get("found"):
                return {
                    "name": res.get("name"),
                    "token": res.get("abbreviation"),
                    "uid": res.get("cas") or uid,
                    "smiles": res.get("smiles"),
                    "inchi_key": res.get("inchi_key"),
                }
        except Exception as e:
            # Fallback: try legacy context lookup if it exists
            try:
                from ..context import lookup
                res = lookup(uid)
                if isinstance(res, dict) and res.get("found") and isinstance(res.get("record"), dict):
                    rec = res["record"]
                    return {
                        "name": rec.get("name"),
                        "token": rec.get("token"),
                        "uid": rec.get("uid") or uid,
                        "smiles": rec.get("smiles"),
                    }
            except Exception:
                pass
        
        # Default fallback: indicate the reagent is unknown
        return {
            "uid": uid,
            "name": f"[{role.replace('_', ' ').title()}] {uid}",  # e.g., "[Base] 7778-53-2"
            "unknown": True
        }
    
    # Helper: Create chemical payload
    def _chemical_payload(uid: Optional[str], role: str) -> Dict[str, Any] | None:
        if not uid:
            return None
        rec = _lookup(uid, role)
        
        # Use name from database, or format unknown reagent nicely
        display_name = rec.get("name")
        if rec.get("unknown"):
            # For unknown reagents, show CAS number with role label
            display_name = f"[Unknown {role.replace('_', ' ')}] CAS {rec.get('uid')}"
        
        return {
            "name": display_name or rec.get("token") or rec.get("uid") or uid,
            "abbreviation": rec.get("token") or None,
            "cas": rec.get("uid") or uid,
            "smiles": rec.get("smiles"),
            "equivalents": None,
            "role": role,
        }
    
    # Reactant chemicals
    reactants_chems: List[Dict[str, Any]] = []
    for r in (norm.get("reactants") or []):
        smi = r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or ""
        reactants_chems.append({
            "name": None,
            "cas": None,
            "smiles": smi or None,
            "equivalents": None,
            "role": "starting_material",
        })
    
    # Catalyst system from precedents
    def _cat_items() -> List[Dict[str, Any]]:
        items: List[Dict[str, Any]] = []
        src = next(
            (p for p in group if p.get("catalytic_system") or p.get("full_system") or p.get("catalyst")),
            None
        )
        fs = None
        if src:
            fs = src.get("catalytic_system")
            if not fs:
                cat = src.get("catalyst") or {}
                fs = src.get("full_system") or (
                    cat.get("full_system") if isinstance(cat, dict) else None
                )
        if not isinstance(fs, list):
            return items
        
        def role_for(name: str) -> str:
            n = (name or "").lower()
            if any(tok in n for tok in ["pd", "palladium", "cu", "copper", "ni", "nickel"]):
                return "metal_precursor"
            return "ligand"
        
        for it in fs:
            nm = (it or {}).get("name")
            cs = (it or {}).get("cas")
            detected_role = role_for(str(nm or ""))
            
            # Use _chemical_payload to get enriched info from database
            if cs:
                # Look up by CAS number
                enriched = _chemical_payload(cs, detected_role)
                if enriched:
                    items.append(enriched)
            elif nm:
                # Look up by name
                enriched = _chemical_payload(nm, detected_role)
                if enriched:
                    items.append(enriched)
            else:
                # No identifier - add placeholder
                items.append({
                    "name": "Unknown catalyst",
                    "cas": None,
                    "smiles": None,
                    "equivalents": None,
                    "role": detected_role,
                })
        return items
    
    cat_items = _cat_items()
    
    # Condition text
    cond_text = {
        "temperature": (f"{int(T_med)} °C" if isinstance(T_med, (int, float)) else None),
        "time": (f"{t_med} h" if isinstance(t_med, (int, float)) else None),
        "atmosphere": None,
    }
    
    core_group_size = len(group) if group else 0
    combo_counts = Counter(
        (str(p.get("base_uid") or ""), str(p.get("solvent_uid") or ""))
        for p in group
    )
    overall_combo_counts = Counter(
        (str(p.get("base_uid") or ""), str(p.get("solvent_uid") or ""))
        for p in precs
    )
    
    # Helper: Build condition variant
    def _build_variant(b_uid: Optional[str], s_uid: Optional[str], rank: int) -> Dict[str, Any]:
        b_key = b_uid or ""
        s_key = s_uid or ""
        
        chems = list(reactants_chems) + list(cat_items)
        base_payload = _chemical_payload(b_uid, "base")
        solvent_payload = _chemical_payload(s_uid, "solvent")
        
        if base_payload:
            chems.append(base_payload)
        if solvent_payload:
            chems.append(solvent_payload)
        
        support_count = combo_counts.get((b_key, s_key), 0)
        if support_count == 0:
            support_count = overall_combo_counts.get((b_key, s_key), 0)
        
        denom = core_group_size if core_group_size else len(precs)
        support_fraction = (support_count / denom) if denom else 0.0
        
        def _matching_precedents(b: str, s: str, src: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
            return [
                p for p in src
                if str(p.get("base_uid") or "") == b and str(p.get("solvent_uid") or "") == s
            ]
        
        matched = _matching_precedents(b_key, s_key, group or precs)
        precedent_examples = [
            {
                "reaction_id": p.get("reaction_id"),
                "reference": p.get("reference"),
                "yield_pct": p.get("yield_pct"),
                "core": p.get("core"),
            }
            for p in matched[:3]
            if p.get("reaction_id")
        ]
        
        summary = {
            "rank": rank,
            "core": chosen_core,
            "base": base_payload,
            "solvent": solvent_payload,
            "confidence": round(float(conf), 3),
            "support": {
                "count": support_count,
                "fraction_core": round(float(support_fraction), 3) if support_fraction else 0.0,
                "reference_population": core_group_size if core_group_size else len(precs),
            },
            "precedents": precedent_examples,
        }
        
        variant = {
            "rank": rank,
            "reaction": {"smiles": norm.get("normalized") or reaction},
            "chemicals": chems,
            "conditions": cond_text,
            "summary": summary,
            "combo": {"base_uid": b_uid, "solvent_uid": s_uid},
        }
        return variant
    
    # Generate variants
    combos: List[Tuple[str, str]] = []
    seen_combos: set[Tuple[str, str]] = set()
    
    def _add_combo(b: Optional[str], s: Optional[str]) -> None:
        key = (b or "", s or "")
        if key in seen_combos:
            return
        seen_combos.add(key)
        combos.append(key)
    
    # Priority: recommended combo first, then by frequency
    _add_combo(base_pick, solv_pick)
    for combo, _ in combo_counts.most_common():
        _add_combo(combo[0], combo[1])
    for combo, _ in overall_combo_counts.most_common():
        _add_combo(combo[0], combo[1])
    for b, _ in base_counts.most_common():
        _add_combo(b, solv_pick)
    for s, _ in solv_counts.most_common():
        _add_combo(base_pick, s)
    
    if not combos:
        _add_combo(None, None)
    
    variants: List[Dict[str, Any]] = []
    for combo in combos[:max_variants]:
        variants.append(_build_variant(combo[0], combo[1], len(variants) + 1))
    
    return {
        "meta": {
            "status": "success",
            "model": "drfp_similarity",
            "version": "2.0",
        },
        "input": {
            "reaction_smiles": norm.get("normalized") or reaction,
            "family": fam,
        },
        "detection": {
            "family": fam,
            "confidence": round(float(conf), 3),
        },
        "recommended_conditions": variants,
        "precedent_summary": {
            "total_precedents": len(precs),
            "chosen_core": chosen_core,
            "core_support": core_group_size,
        },
        "precedents_used": _build_precedent_details(precs, chosen_core, group),
    }
