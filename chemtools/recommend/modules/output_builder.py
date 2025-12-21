"""
Formatted output builder for recommendation results.

This module builds API-friendly multi-variant outputs with detailed chemical
and precedent information.
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional, Tuple
from collections import Counter

from ..utils import friendly_family_label
from .precedent_builder import build_precedent_details


def build_formatted_output(
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
    
    Args:
        norm: Normalized reaction dict from normalize_reaction()
        reaction: Original reaction SMILES
        fam: Detected reaction family
        recommendation: Primary recommendation dict
        alternatives: Alternative cores/bases/solvents
        precs: All precedents from search
        group: Precedents matching chosen core
        chosen_core: Chosen catalyst core
        base_pick: Chosen base
        solv_pick: Chosen solvent
        T_med: Median temperature
        t_med: Median time
        conf: Confidence score
        base_counts: Counter of base frequencies
        solv_counts: Counter of solvent frequencies
        max_variants: Maximum number of variants to generate
        
    Returns:
        Formatted output dict with multi-variant recommendations
    """
    # Helper: Lookup reagent info from reagent database
    def _lookup(uid: str, role: str) -> Dict[str, Any]:
        """Look up reagent information by CAS number or name."""
        try:
            from ... import reagent
            
            # Map role to reagent type for lookup
            role_to_type = {
                "base": "base",
                "solvent": "solvent",
                "ligand": "ligand",
                "metal_catalyst": "metal_catalyst",
                "catalyst": "metal_catalyst",
            }
            reagent_type = role_to_type.get(role, role)
            
            # Try to enrich the reagent info
            res = reagent.enrich_reagent_info(uid, reagent_type)
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
                from ...context import lookup
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
        """
        Extract catalyst system from precedents.
        
        Strategy: Use the catalytic_system from the top-ranked precedent overall,
        preferring specific catalyst complexes over generic metals. This preserves
        actionable catalyst information (e.g., "Pd(dppf)Cl2" instead of generic "Pd").
        
        Falls back to group-specific precedents if top precedent has no catalytic_system.
        """
        items: List[Dict[str, Any]] = []
        
        # Strategy: Try top-ranked precedent first (best similarity match)
        # This often has the most specific catalyst information
        src = None
        
        # 1. Try top-ranked precedent overall (highest similarity)
        if precs:
            top_prec = precs[0]
            if top_prec.get("catalytic_system") or top_prec.get("full_system") or top_prec.get("catalyst"):
                src = top_prec
        
        # 2. Fall back to first precedent in chosen core group
        if not src and group:
            src = next(
                (p for p in group if p.get("catalytic_system") or p.get("full_system") or p.get("catalyst")),
                None
            )
        
        # 3. Fall back to any precedent with catalytic_system
        if not src and precs:
            src = next(
                (p for p in precs if p.get("catalytic_system") or p.get("full_system") or p.get("catalyst")),
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
            if any(tok in n for tok in ["pd", "palladium", "cu", "copper", "ni", "nickel", "ir", "iridium", "ru", "ruthenium"]):
                return "metal_catalyst"
            return "ligand"
        
        for it in fs:
            nm = (it or {}).get("name")
            cs = (it or {}).get("cas")
            detected_role = role_for(str(nm or ""))
            
            # Try database lookup for additional info (SMILES, etc.)
            enriched = None
            if cs:
                enriched = _chemical_payload(cs, detected_role)
            elif nm:
                enriched = _chemical_payload(nm, detected_role)
            
            if enriched:
                # If lookup marked as "unknown", override with precedent's original name
                if enriched.get("unknown") or "[Unknown" in str(enriched.get("name", "")):
                    if nm:
                        # Use the specific catalyst name from the precedent
                        enriched["name"] = nm
                        enriched.pop("unknown", None)  # Remove unknown flag
                items.append(enriched)
            elif nm or cs:
                # No database match - use original information from precedent
                items.append({
                    "name": nm or f"CAS {cs}",
                    "cas": cs,
                    "smiles": None,
                    "equivalents": None,
                    "role": detected_role,
                })
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
    
    # Fallback: If no catalytic_system found but we have a chosen_core,
    # add a generic metal entry to ensure recommendations include at least the core metal
    if not cat_items and chosen_core:
        # Try to create a generic metal entry from the core string
        # e.g., "Pd" -> lookup "Pd" or "Palladium"
        core_payload = _chemical_payload(chosen_core, "metal_catalyst")
        if core_payload:
            cat_items.append(core_payload)
    
    # Condition text
    cond_text = {
        "temperature": (f"{int(T_med)} °C" if isinstance(T_med, (int, float)) else None),
        "time": (f"{t_med} h" if isinstance(t_med, (int, float)) else None),
        "atmosphere": None,
    }
    
    # Add cross-family metadata if available (from first precedent in group or all precs)
    source_precs = group if group else precs
    if source_precs:
        first_prec = source_precs[0]
        if first_prec.get("cross_family_metadata"):
            cond_text["cross_family_metadata"] = first_prec.get("cross_family_metadata")
        # Also add individual compatibility fields for easier access
        if first_prec.get("mechanism_similarity") is not None:
            cond_text["mechanism_similarity"] = first_prec.get("mechanism_similarity")
        if first_prec.get("mechanism_status"):
            cond_text["mechanism_status"] = first_prec.get("mechanism_status")
        if first_prec.get("reaction_type_status"):
            cond_text["reaction_type_status"] = first_prec.get("reaction_type_status")
        if first_prec.get("rxn_type"):
            cond_text["reaction_family"] = first_prec.get("rxn_type")
    
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
                "rxn_type": p.get("rxn_type"),  # Include reaction family/type for dataset identification
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
            "family_label": friendly_family_label(fam),
        },
        "detection": {
            "family": fam,
            "family_label": friendly_family_label(fam),
            "confidence": round(float(conf), 3),
        },
        "recommended_conditions": variants,
        "precedent_summary": {
            "total_precedents": len(precs),
            "chosen_core": chosen_core,
            "core_support": core_group_size,
        },
        "precedents_used": build_precedent_details(precs, chosen_core, group),
    }
