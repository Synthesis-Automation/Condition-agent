"""
Main DRFP-based recommendation engine.

This module contains the primary recommendation function that uses reaction
fingerprinting (DRFP) for similarity-based precedent search.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional
from collections import Counter

from ...smiles import normalize_reaction
from ...reaction_type_detection import detect_reaction  # New unified API
from ...router import resolve_reaction_family
from ... import precedent, explain
from ...taxonomy import reaction_catalog as _reaction_catalog
from ..utils import canonical_family, median, pick_with_constraints

# No longer needed - use unified detect_reaction() API instead
# Old rxn-insight integration removed

# New analysis module integration for improved reaction type detection and reactant classification
try:
    from ...featurizers.analysis.reaction_context import classify_reactants_with_context
    _HAS_ANALYSIS_MODULE = True
except Exception:
    _HAS_ANALYSIS_MODULE = False


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
    search_all_families: bool = False,
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
        family_override: Override auto-detected reaction family (ignored if search_all_families=True)
        max_variants: Maximum number of condition variants to generate (default: 3)
        rerank_strategy: Precedent reranking strategy (default: 'rule')
            - 'rule': Boost precedents matching chemical rules (e.g., Ullmann→Cu)
            - 'analytics': Boost precedents using popular reagents from dataset
            - 'none': Use DRFP similarity only (no reranking)
        filter_unknown_reagents: If True, remove precedents with unknown base/solvent reagents (default: False)
            Note: Only filters base and solvent, not catalyst cores (cores are complex "Metal/Ligand" strings)
        search_all_families: If True, search across all reaction families (ignore reaction type) (default: False)
            Note: When True, family detection is skipped and precedents from all families are considered
    
    Returns:
        Dict with keys:
            - input: Normalized reaction SMILES
            - family: Detected reaction family (or "All" if search_all_families=True)
            - features: Feature dict (usually empty when DRFP is used)
            - recommendation: Primary recommendation (core, base, solvent, T_C, time_h, confidence)
            - alternatives: Alternative cores, bases, solvents
            - precedent_pack: Full precedent search results
            - reasons: Human-readable explanation
            - formatted: Formatted output with multiple variants
    
    Example:
        >>> from chemtools.recommend import recommend_from_reaction
        >>> # Standard family-specific search
        >>> results = recommend_from_reaction(
        ...     reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        ...     k=50,
        ...     rerank_strategy='rule',
        ...     filter_unknown_reagents=True
        ... )
        >>> print(results['recommendation']['core'])
        'Cu/phen'  # Correctly identifies Cu for Ullmann, not Pd
        >>> 
        >>> # Cross-family search (search all datasets)
        >>> results = recommend_from_reaction(
        ...     reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        ...     k=100,
        ...     search_all_families=True
        ... )
        >>> print(results['family'])
        'All'  # Searches across all reaction types
    """
    relax = dict(relax or {})
    max_variants = max(1, int(max_variants or 1))
    fam_override_clean = (family_override.strip() if family_override else None)
    if fam_override_clean:
        resolved_override = resolve_reaction_family(fam_override_clean)
        if resolved_override:
            fam_override_clean = resolved_override

    # 1) Normalize reaction SMILES
    norm = normalize_reaction(reaction)
    reactants = [
        (r.get("smiles_norm") or r.get("largest_smiles") or r.get("input") or "")
        for r in (norm.get("reactants") or [])
    ]

    # 2) Detect reaction family
    rxn_smiles_norm = norm.get("normalized") or reaction
    
    # Always try to detect the reaction type (even for cross-family search)
    # Detection helps with filtering and metadata, even when searching all families
    detection_source = "user_supplied" if fam_override_clean else "auto"
    fam = fam_override_clean or "Unknown"
    
    # Priority 1: New analysis module (most comprehensive - includes reactant classification)
    auto_family = None
    analysis_result = None
    use_analysis_module = relax.get("use_analysis_module", True)  # Enabled by default
    
    if use_analysis_module and _HAS_ANALYSIS_MODULE and not fam_override_clean:
        try:
            # Use Two-Pass Approach with auto-detection
            analysis_result = classify_reactants_with_context(
                reaction_smiles=rxn_smiles_norm,
                reaction_type=None,  # Auto-detect
                auto_detect=True
            )
            if analysis_result and analysis_result.reaction_type:
                auto_family = analysis_result.reaction_type
                # Don't store in relax (needs to be hashable for caching)
                # Will be stored separately in the result
        except Exception as e:
            # Fallback to other methods if analysis module fails
            import warnings
            warnings.warn(
                f"Analysis module failed: {e}. Falling back to rxn-insight/rule-based detection.",
                UserWarning
            )
    
    # Priority 2: rxn-insight ML-based detection (if analysis module didn't work)
    use_rxn_insight = relax.get("use_rxn_insight")
    if use_rxn_insight is None:
        import os as _os
        env_off = (_os.environ.get("CHEMTOOLS_DISABLE_RXN_INSIGHT", "").strip().lower() 
                   in {"1", "true", "yes", "on"})
        use_rxn_insight = not env_off
    
    # ML detection removed - use unified detect_reaction() API
    # Old code using _rxn_detect has been removed
    
    # Priority 3: Rule-based detection as final fallback
    # Convert reactants to pseudo-reaction for unified API
    pseudo_reaction = ".".join(reactants) + ">>"
    detection_result = detect_reaction(pseudo_reaction, use_ml=False)
    rule_family = detection_result.get("family") or "Unknown"
    
    # Prioritize: user override > analysis module > rule-based
    if not fam_override_clean:
        detected_id = auto_family or rule_family or "Unknown"
    else:
        detected_id = fam_override_clean
    fam = detected_id
    
    # Normalize to canonical precedent database family name
    # (e.g., ullmann_cn -> C_N_Coupling, suzuki_miyaura -> Suzuki)
    canonical_fam = canonical_family(fam)
    
    # For cross-family search: detect but search all families
    if search_all_families:
        # Keep detection metadata but set fam to None for cross-family precedent search
        detection_source = "cross_family_search"
        detected_fam = canonical_fam  # Store dataset family for metadata
        fam = None  # Signal for cross-family search in precedent.knn()
    else:
        detected_fam = canonical_fam
        fam = canonical_fam
    
    # 3) Features: Include taxonomy-aware metadata for similarity fallback
    features: Dict[str, Any] = {}
    if detected_id and detected_id != "Unknown":
        features["reaction_type"] = detected_id
        resolved_id = _reaction_catalog.resolve_reaction_type(detected_id)
        if resolved_id is None:
            resolved_id = _reaction_catalog.resolve_reaction_type(str(detected_id).lower())
        if resolved_id:
            record = _reaction_catalog.get_reaction_type(resolved_id)
            if record and record.category:
                features["reaction_category"] = record.category
    
    # 4) Retrieve precedents using DRFP-based k-NN search
    relax.setdefault("reaction_smiles", rxn_smiles_norm)
    
    # For cross-family search, check if unified DRFP index exists
    if search_all_families and "use_drfp" not in relax:
        # Check if unified index exists
        try:
            from ...util.drfp_storage import get_unified_drfp_path
            import os
            unified_path = get_unified_drfp_path()
            if os.path.exists(unified_path):
                # Unified index available - enable DRFP!
                relax.setdefault("use_drfp", True)
                if relax.get("debug_timing", False):
                    import warnings
                    warnings.warn(
                        f"Cross-family search with DRFP enabled (unified index found at {unified_path})",
                        UserWarning
                    )
            else:
                # No unified index - disable DRFP to avoid slow on-demand computation
                relax.setdefault("use_drfp", False)
                import warnings
                warnings.warn(
                    "Cross-family search with DRFP disabled (no unified index found). "
                    "Results will use feature-based similarity only. "
                    "To enable DRFP, run: python scripts/build_unified_drfp_index.py",
                    UserWarning
                )
        except ImportError:
            # DRFP storage utilities not available
            relax.setdefault("use_drfp", False)
            import warnings
            warnings.warn(
                "Cross-family search with DRFP disabled (no precomputed fingerprints available). "
                "Results will use feature-based similarity only.",
                UserWarning
            )
    else:
        relax.setdefault("use_drfp", True)  # Enable DRFP by default for family-specific search
    
    relax.setdefault("precompute_drfp", False)  # Use binary NPZ files (faster)
    relax.setdefault("selective_loading", True)  # Load only needed family
    relax.setdefault("filter_by_reagent_database", bool(filter_unknown_reagents))
    
    pack = precedent.knn(family=fam, features=features, k=int(k), relax=relax)
    
    precs: List[Dict[str, Any]] = list(pack.get("precedents") or [])
    support = int(pack.get("support") or len(precs))
    
    # Initialize reasoning list for tracking filtering and reranking
    rerank_reasons = []
    
    # 5) Optional filtering: Remove precedents with unknown reagents
    if filter_unknown_reagents and precs and not relax.get("filter_by_reagent_database", False):
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
    # Skip rule-based reranking for cross-family search (family-specific rules don't apply)
    if rerank_strategy in ['rule', 'analytics'] and precs and not search_all_families:
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
    elif search_all_families and rerank_strategy != 'none':
        rerank_reasons.append(f"Reranking skipped for cross-family search ('{rerank_strategy}' not applicable)")
    
    # 6.5) Cross-family mechanism-aware marking and reranking
    if search_all_families and precs and detected_fam:  # Use detected_fam instead of fam
        try:
            from ...precedent.mechanism_similarity import enhance_precedent_similarity, calculate_mechanism_similarity
            
            # Get configuration parameters
            reaction_type_threshold = float(relax.get('reaction_type_threshold', 0.15))
            mechanism_threshold = float(relax.get('mechanism_similarity_threshold', 0.4))
            mechanism_weight = float(relax.get('mechanism_weight', 0.3))
            
            # Step 1: Mark reaction type compatibility (don't filter, just mark)
            if reaction_type_threshold > 0:
                from ...precedent.mechanism_similarity import get_representation_status
                
                type_counts = Counter([p.get('rxn_type', 'unknown') for p in precs])
                total_precs = len(precs)
                
                marked_count = 0
                for prec in precs:
                    prec_type = prec.get('rxn_type', 'unknown')
                    type_fraction = type_counts[prec_type] / total_precs if total_precs > 0 else 0
                    
                    # Mark all precedents with their reaction type compatibility
                    prec['reaction_type_fraction'] = type_fraction
                    prec['reaction_type_status'] = get_representation_status(type_fraction, reaction_type_threshold)
                    
                    if prec['reaction_type_status'] == 'underrepresented':
                        marked_count += 1
                
                if marked_count > 0:
                    rerank_reasons.append(f"Marked {marked_count} precedents as underrepresented reaction types (threshold: {reaction_type_threshold:.1%})")
            
            # Step 2: Calculate mechanism similarity and mark compatibility
            from ...precedent.mechanism_similarity import get_compatibility_status
            
            mechanism_marked = 0
            for prec in precs:
                prec_family = prec.get('rxn_type', '')
                mechanism_sim = calculate_mechanism_similarity(detected_fam, prec_family)  # Use detected_fam
                prec['mechanism_similarity'] = mechanism_sim
                prec['mechanism_status'] = get_compatibility_status(mechanism_sim)
                
                if prec['mechanism_status'] != 'compatible':
                    mechanism_marked += 1
            
            if mechanism_marked > 0:
                rerank_reasons.append(f"Marked {mechanism_marked} precedents with low/moderate mechanism compatibility")
            
            # Step 3: Enhance similarity scores with mechanism awareness and compatibility penalties
            if mechanism_weight > 0 and precs:
                for i, prec in enumerate(precs):
                    base_sim = prec.get('drfp_similarity', prec.get('similarity', 0.5))
                    
                    # Apply mechanism enhancement
                    enhanced_sim = enhance_precedent_similarity(
                        prec, detected_fam, base_sim, mechanism_weight  # Use detected_fam
                    )
                    
                    # Apply compatibility penalties for ranking
                    final_sim = enhanced_sim
                    
                    # Penalty for underrepresented reaction types
                    if prec.get('reaction_type_status') == 'underrepresented':
                        final_sim *= 0.85  # 15% penalty
                    
                    # Penalty for incompatible mechanisms
                    if prec.get('mechanism_status') == 'incompatible':
                        final_sim *= 0.70  # 30% penalty
                    elif prec.get('mechanism_status') == 'moderate':
                        final_sim *= 0.90  # 10% penalty
                    
                    prec['similarity'] = final_sim
                    prec['base_similarity'] = base_sim  # Keep original for reference
                    prec['enhanced_similarity'] = enhanced_sim  # Keep enhanced for reference
                    prec['mechanism_enhanced'] = True
                
                # Re-sort by final similarity (with penalties applied)
                precs.sort(key=lambda p: p.get('similarity', 0), reverse=True)
                rerank_reasons.append(f"Enhanced similarity with mechanism awareness and compatibility penalties (weight: {mechanism_weight:.1f})")
            
            # Step 4: Add cross-family metadata for transparency
            for prec in precs:
                prec['cross_family_metadata'] = {
                    'detected_family': detected_fam,  # Use detected_fam
                    'precedent_family': prec.get('rxn_type', 'unknown'),
                    'mechanism_similarity': prec.get('mechanism_similarity', 0),
                    'mechanism_status': prec.get('mechanism_status', 'unknown'),
                    'reaction_type_status': prec.get('reaction_type_status', 'unknown'),
                    'reaction_type_fraction': prec.get('reaction_type_fraction', 0)
                }
            
        except Exception as e:
            import warnings
            warnings.warn(f"Cross-family mechanism marking failed: {e}. Using all precedents without marking.")
            rerank_reasons.append(f"Mechanism marking error: {e} - using all precedents without marking")
    
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
        "detected_family": detected_fam,  # Dataset family used for precedent search
        "detected_family_id": detected_id,  # Taxonomy family ID (if detected)
        "source": detection_source,
        "search_mode": "cross_family" if search_all_families else "family_specific",
        "analysis_module_used": analysis_result is not None,
        "reactant_classification": None,  # Will be populated if analysis module was used
    }
    
    # Add reactant classification info from analysis module
    if analysis_result is not None:
        try:
            from ...featurizers.analysis.reaction_context import get_reactant_summary
            detection_meta["reactant_classification"] = get_reactant_summary(analysis_result)
        except Exception:
            pass  # Ignore serialization errors
    
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
        "family": "All" if search_all_families else fam,  # Use "All" for cross-family, canonical name for family-specific
        "detected_family": detected_fam,  # The actual detected family (for display/metadata)
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
        "search_all_families": search_all_families,
    }
