"""
Fusion recommender format adapter.

This module converts fusion recommender output (multi-source evidence)
to the core.py format for API compatibility.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional
from collections import Counter


def convert_fusion_to_core_format(
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
        
    Args:
        fusion_results: Output from fusion recommender
        reaction: Original reaction SMILES
        rxn_smiles_norm: Normalized reaction SMILES
        fam: Detected reaction family
        features: Feature dict
        detection_meta: Detection metadata
        
    Returns:
        Dict in core.py format with all expected keys
    """
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
    
    # Build alternatives from all recommendations (using converted dicts)
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
    
    formatted = build_formatted_output_from_fusion(
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


def build_formatted_output_from_fusion(
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
    
    Args:
        fusion_results: Converted fusion results (dicts, not dataclasses)
        reaction: Original reaction SMILES
        rxn_smiles_norm: Normalized reaction SMILES
        fam: Detected reaction family
        recommendation: Primary recommendation
        alternatives: Alternative reagents
        precedents: Precedent reactions
        weights: Adaptive weights from fusion
        
    Returns:
        Formatted output dict matching standard API format
    """
    from ...smiles import normalize_reaction
    
    norm = normalize_reaction(reaction)
    recommended_conditions = fusion_results.get('recommended_conditions', [])
    
    # Helper: Lookup reagent info
    def _lookup(uid: str, role: str) -> Dict[str, Any]:
        """Look up reagent information by CAS number or name."""
        try:
            from ... import reagent
            
            role_to_type = {
                "base": "base",
                "solvent": "solvent",
                "ligand": "ligand",
                "metal_catalyst": "metal_catalyst",
                "catalyst": "metal_catalyst",
            }
            reagent_type = role_to_type.get(role, role)
            
            res = reagent.enrich_reagent_info(uid, reagent_type)
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
