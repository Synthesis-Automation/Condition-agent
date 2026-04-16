"""
Composite facade tools for ChemCoworker.

These tools orchestrate multiple low-level tools into higher-level operations so
the LLM can solve common tasks with fewer tool calls and less schema clutter.
Internal tools remain available for debugging and power-user workflows.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from ._base import ToolPlugin
from ._helpers import _error, _success, _validate_reaction_smiles


_CONDITION_SOURCE_MODES = ("literature", "motif", "similarity", "rules")
_CONDITION_SOURCE_WEIGHTS: Dict[str, float] = {
    "literature": 1.00,
    "motif": 0.95,
    "similarity": 0.75,
    "rules": 0.45,
}


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None:
            return default
        return float(value)
    except Exception:
        return default


def _as_int(value: Any, default: int = 0) -> int:
    try:
        if value is None:
            return default
        return int(value)
    except Exception:
        return default


def _condition_fingerprint(rec: Dict[str, Any]) -> str:
    """Stable fingerprint for deduplicating condition sets across sources."""
    fields = [
        "catalyst",
        "ligand",
        "base",
        "solvent",
        "secondary_solvent",
        "additive",
        "coupling_reagent",
        "temperature",
        "atmosphere",
    ]
    parts = []
    for key in fields:
        val = str(rec.get(key) or "").strip().lower()
        parts.append(f"{key}={val}")
    return "|".join(parts)


def _entry_support_score(rec: Dict[str, Any]) -> float:
    """Deterministic quality score from support metrics (0..~1.2)."""
    confidence = max(0.0, min(1.0, _as_float(rec.get("confidence"), 0.0)))
    success_rate = max(0.0, min(1.0, _as_float(rec.get("success_rate"), 0.0)))
    avg_yield = _as_float(rec.get("avg_yield"), 0.0)
    if avg_yield > 1.0:
        avg_yield = avg_yield / 100.0
    avg_yield = max(0.0, min(1.0, avg_yield))
    median_yield = _as_float(rec.get("median_yield"), 0.0)
    if median_yield > 1.0:
        median_yield = median_yield / 100.0
    median_yield = max(0.0, min(1.0, median_yield))
    n_exp = _as_int(rec.get("num_experiments"), 0)
    n_exp_norm = max(0.0, min(1.0, n_exp / 20.0))
    return round(
        0.35 * confidence
        + 0.20 * success_rate
        + 0.20 * avg_yield
        + 0.10 * median_yield
        + 0.15 * n_exp_norm,
        4,
    )


def _recommend_conditions_with_strategy(
    reaction_smiles: str,
    top_k: int = 5,
    condition_strategy: str = "auto",
    condition_source_mode: str = "",
    condition_selection_mode: str = "best",
) -> Dict[str, Any]:
    """Run recommend_conditions in single-source or full multi-source mode."""
    from .conditions import _compose_condition_candidates, _recommend_conditions

    reaction_smiles, rxn_err = _validate_reaction_smiles(reaction_smiles, require_product=True)
    if rxn_err:
        return _error(f"Invalid reaction_smiles: {rxn_err}")

    top_k = max(1, min(int(top_k or 5), 10))
    strategy = (condition_strategy or "auto").strip().lower()
    source_mode = (condition_source_mode or "").strip().lower()
    if source_mode == "all":
        source_mode = ""

    if strategy not in {"auto", "single", "full", "balanced"}:
        return _error("condition_strategy must be one of: auto, single, full")
    if strategy == "balanced":
        strategy = "full"
    if source_mode and source_mode not in _CONDITION_SOURCE_MODES:
        return _error("condition_source_mode must be one of: literature, motif, similarity, rules, all")

    # Auto defaults to a single recommend_conditions call; callers can opt into
    # cross-source full ensemble mode explicitly via condition_strategy='full'.
    if strategy == "auto":
        strategy = "single"

    if strategy == "single":
        result = _recommend_conditions(
            reaction_smiles=reaction_smiles,
            top_k=top_k,
            source_group=source_mode or "",
            selection_mode=condition_selection_mode,
        )
        if isinstance(result, dict):
            out = dict(result)
            out["condition_strategy"] = "single"
            out["condition_source_mode"] = source_mode or "all"
            return out
        return _error("Unexpected recommend_conditions result")

    composed = _compose_condition_candidates(
        reaction_smiles=reaction_smiles,
        top_k=top_k,
        sources=list(_CONDITION_SOURCE_MODES),
        selection_mode=condition_selection_mode,
    )
    if isinstance(composed, dict):
        out = dict(composed)
        out["condition_strategy"] = "full"
        out["condition_source_mode"] = "full"
        out["source_group"] = "full"
        return out
    return _error("Unexpected compose_condition_candidates result")


def _resolve_chemical(
    identifier: str = "",
    smiles: str = "",
    mode: str = "auto",
) -> Dict[str, Any]:
    """Unified name/CAS/InChI<->SMILES helper.

    Args:
        identifier: Name, CAS, or InChI to resolve to SMILES.
        smiles: SMILES to annotate with names/CAS metadata.
        mode: 'auto', 'to_smiles', or 'from_smiles'.
    """
    from .name_resolver import _resolve_to_smiles, _smiles_to_info

    mode_norm = (mode or "auto").strip().lower()
    if mode_norm not in {"auto", "to_smiles", "from_smiles"}:
        return _error("mode must be one of: auto, to_smiles, from_smiles")

    if mode_norm == "auto":
        if identifier.strip() and not smiles.strip():
            mode_norm = "to_smiles"
        elif smiles.strip() and not identifier.strip():
            mode_norm = "from_smiles"
        elif identifier.strip() and smiles.strip():
            mode_norm = "to_smiles"
        else:
            return _error("Provide either identifier or smiles")

    result = (
        _resolve_to_smiles(identifier=identifier)
        if mode_norm == "to_smiles"
        else _smiles_to_info(smiles=smiles)
    )
    if isinstance(result, dict):
        out = dict(result)
        out["mode_used"] = mode_norm
        out["tool_used"] = "resolve_to_smiles" if mode_norm == "to_smiles" else "smiles_to_info"
        return out
    return _error("Unexpected resolver result")


def _reagent_assistant(
    mode: str = "auto",
    query: str = "",
    role: str = "",
    family: str = "",
) -> Dict[str, Any]:
    """Unified reagent lookup/list facade."""
    from .reagent import _list_reagents_by_role, _lookup_reagent

    mode_norm = (mode or "auto").strip().lower()
    if mode_norm not in {"auto", "lookup", "list"}:
        return _error("mode must be one of: auto, lookup, list")

    if mode_norm == "auto":
        mode_norm = "list" if role.strip() else "lookup"

    if mode_norm == "lookup":
        if not query.strip():
            return _error("query is required for mode='lookup'")
        result = _lookup_reagent(query=query)
        tool_used = "lookup_reagent"
    else:
        if not role.strip():
            return _error("role is required for mode='list'")
        result = _list_reagents_by_role(role=role, family=family or None)
        tool_used = "list_reagents_by_role"

    if isinstance(result, dict):
        out = dict(result)
        out["mode_used"] = mode_norm
        out["tool_used"] = tool_used
        return out
    return _error("Unexpected reagent result")


def _search_knowledge(
    query: str,
    top_k: int = 5,
    include_notes: bool = True,
    include_literature: bool = True,
    note_types: Optional[List[str]] = None,
    reaction_note_id: str = "",
) -> Dict[str, Any]:
    """Search notes and literature together, with optional targeted reaction note."""
    from .literature import _search_literature
    from .notes import _read_notes, _search_notes

    if not query.strip() and not reaction_note_id.strip():
        return _error("Provide query and/or reaction_note_id")

    top_k = max(1, min(int(top_k or 5), 10))
    per_source_top_k = max(1, min(5, top_k))

    notes_result: Dict[str, Any] = {"success": True, "found": 0, "results": []}
    literature_result: Dict[str, Any] = {"success": True, "found": 0, "results": []}
    reaction_note: Dict[str, Any] = {}

    if include_notes and query.strip():
        notes_result = _search_notes(
            query=query,
            note_types=note_types,
            top_k=per_source_top_k,
        )
    if include_literature and query.strip():
        literature_result = _search_literature(query=query, top_k=per_source_top_k)
    if reaction_note_id.strip():
        note = _read_notes(id=reaction_note_id.strip(), note_type="reactions", max_chars=2500)
        if note.get("success") and note.get("has_notes"):
            reaction_note = {
                "id": note.get("id"),
                "note_type": note.get("note_type", "reactions"),
                "char_count": note.get("char_count", 0),
                "content_preview": (note.get("content") or "")[:1500],
                "notes_file": note.get("notes_file", ""),
            }
        elif note.get("success"):
            reaction_note = {
                "id": note.get("id", reaction_note_id),
                "has_notes": False,
                "note": note.get("note", "No notes found"),
            }

    notes_found = int(notes_result.get("found", 0) or 0) if notes_result.get("success", True) else 0
    lit_found = int(literature_result.get("found", 0) or 0) if literature_result.get("success", True) else 0
    total_found = notes_found + lit_found + (1 if reaction_note else 0)

    return _success({
        "query": query,
        "top_k": top_k,
        "sources_searched": {
            "notes": bool(include_notes),
            "literature": bool(include_literature),
        },
        "source_counts": {
            "notes_hits": notes_found,
            "literature_hits": lit_found,
            "targeted_reaction_note": 1 if reaction_note else 0,
        },
        "total_found": total_found,
        "notes": notes_result,
        "literature": literature_result,
        "reaction_note": reaction_note,
    })


def _read_knowledge(
    id_or_filename: str,
    source: str = "auto",
    note_type: str = "",
    max_chars: int = 8000,
    fetch_if_url: bool = True,
) -> Dict[str, Any]:
    """Unified note/literature reader, optionally fetching a URL directly."""
    from .literature import _fetch_webpage, _read_literature_source
    from .notes import _read_notes

    target = (id_or_filename or "").strip()
    if not target:
        return _error("id_or_filename is required")

    source_norm = (source or "auto").strip().lower()
    if source_norm not in {"auto", "notes", "literature", "url"}:
        return _error("source must be one of: auto, notes, literature, url")

    is_url = target.startswith(("http://", "https://"))
    if source_norm == "auto":
        if is_url and fetch_if_url:
            source_norm = "url"
        elif target.lower().endswith((".txt", ".md", ".pdf")) or "/" in target or "\\" in target:
            source_norm = "literature"
        else:
            source_norm = "notes"

    if source_norm == "url":
        result = _fetch_webpage(url=target, max_chars=max_chars)
        tool_used = "fetch_webpage"
    elif source_norm == "literature":
        result = _read_literature_source(filename=target, max_chars=max_chars)
        tool_used = "read_literature_source"
    else:
        result = _read_notes(id=target, note_type=(note_type or None), max_chars=max_chars)
        tool_used = "read_notes"

    if isinstance(result, dict):
        out = dict(result)
        out["source_used"] = source_norm
        out["tool_used"] = tool_used
        return out
    return _error("Unexpected knowledge read result")


def _analyze_reaction(
    reaction_smiles: str,
    include_conditions: bool = False,
    conditions_top_k: int = 5,
    include_notes: bool = False,
    include_bond_changes: bool = True,
    condition_strategy: str = "auto",
    condition_source_mode: str = "",
    condition_selection_mode: str = "best",
) -> Dict[str, Any]:
    """High-level reaction analysis: normalize, classify, inspect, and optionally rank conditions."""
    from .chemistry import (
        _analyze_bond_changes,
        _detect_reaction_type,
        _inspect_functional_groups,
        _normalize_reaction,
    )
    from .notes import _read_notes

    norm = _normalize_reaction(reaction_smiles)
    if not norm.get("success"):
        return norm
    if not norm.get("is_reaction", False):
        return _error("analyze_reaction expects a reaction SMILES in 'reactants>>products' format")

    det = _detect_reaction_type(reaction_smiles=reaction_smiles)
    bond = _analyze_bond_changes(reaction_smiles=reaction_smiles) if include_bond_changes else {}

    reactant_profiles: List[Dict[str, Any]] = []
    product_profiles: List[Dict[str, Any]] = []
    for smi in (norm.get("reactants") or [])[:2]:
        fg = _inspect_functional_groups(smiles=smi)
        reactant_profiles.append({
            "smiles": smi,
            "detected_groups": fg.get("detected_groups", [])[:20] if fg.get("success") else [],
            "categories": fg.get("categories", {}) if fg.get("success") else {},
        })
    for smi in (norm.get("products") or [])[:2]:
        fg = _inspect_functional_groups(smiles=smi)
        product_profiles.append({
            "smiles": smi,
            "detected_groups": fg.get("detected_groups", [])[:20] if fg.get("success") else [],
            "categories": fg.get("categories", {}) if fg.get("success") else {},
        })

    conditions_result: Dict[str, Any] = {}
    if include_conditions:
        conditions_result = _recommend_conditions_with_strategy(
            reaction_smiles=reaction_smiles,
            top_k=max(1, min(int(conditions_top_k or 5), 10)),
            condition_strategy=condition_strategy,
            condition_source_mode=condition_source_mode,
            condition_selection_mode=condition_selection_mode,
        )

    reaction_note: Dict[str, Any] = {}
    if include_notes and det.get("success"):
        rxn_id = str(det.get("reaction_type_id") or det.get("reaction_type") or "").strip()
        if rxn_id and rxn_id.lower() != "unknown":
            note = _read_notes(id=rxn_id, note_type="reactions", max_chars=2000)
            if note.get("success") and note.get("has_notes"):
                reaction_note = {
                    "id": note.get("id"),
                    "char_count": note.get("char_count", 0),
                    "content_preview": (note.get("content") or "")[:1200],
                }

    return _success({
        "reaction_smiles": reaction_smiles,
        "normalization": norm,
        "reaction_type": det,
        "bond_changes": bond,
        "reactant_functional_groups": reactant_profiles,
        "product_functional_groups": product_profiles,
        "conditions": conditions_result,
        "reaction_note": reaction_note,
        "summary": {
            "reaction_type_id": det.get("reaction_type_id") if isinstance(det, dict) else "",
            "reaction_type_confidence": det.get("confidence") if isinstance(det, dict) else None,
            "top_condition_count": len((conditions_result or {}).get("recommendations", []) or []),
        },
    })


def _recommend_reaction_conditions(
    reaction_smiles: str,
    top_k: int = 5,
    condition_strategy: str = "auto",
    condition_source_mode: str = "",
    condition_selection_mode: str = "best",
) -> Dict[str, Any]:
    """Explicit facade for reaction-condition recommendation."""
    rxn, rxn_err = _validate_reaction_smiles(reaction_smiles, require_product=True)
    if rxn_err:
        return _error(f"Invalid reaction_smiles: {rxn_err}")
    return _recommend_conditions_with_strategy(
        reaction_smiles=rxn,
        top_k=max(1, min(int(top_k or 5), 10)),
        condition_strategy=condition_strategy,
        condition_source_mode=condition_source_mode,
        condition_selection_mode=condition_selection_mode,
    )


def _retrosynthesis_step(
    target_smiles: str,
    top_k: int = 3,
    include_precedent: bool = False,
    include_conditions: bool = True,
    validate_top_n: int = 1,
    condition_strategy: str = "auto",
    condition_source_mode: str = "",
    condition_selection_mode: str = "best",
) -> Dict[str, Any]:
    """Single-call retrosynthesis pass with automatic fallback cascade.

    Pipeline: inspect → retrons → disconnections → (fallback: HTE templates,
    product similarity) → validation → conditions → precedent.

    When standard retrons yield 0 disconnections the cascade fires:
      1. apply_hte_templates (35+ atom-mapped SMARTS)
      2. search_by_product_similarity (Morgan FP against 231k HTE reactions)
    The best data-driven disconnection from templates or product-similarity
    is promoted to `top_disconnection` so conditions are always grounded.
    """
    from .reaction_eval import _check_retro_consistency
    from .retrosynthesis import (
        _apply_hte_templates,
        _find_retro_precedent,
        _generate_disconnections,
        _identify_retrons,
        _inspect_target,
        _search_by_product_similarity,
    )

    top_k = max(1, min(int(top_k or 3), 10))
    validate_top_n = max(0, min(int(validate_top_n or 0), 3))

    def _has_concrete_precursors(row: Dict[str, Any]) -> bool:
        return bool(str(row.get("precursor_1") or "").strip() and str(row.get("precursor_2") or "").strip())

    def _with_reaction_smiles(row: Dict[str, Any]) -> Dict[str, Any]:
        out = dict(row or {})
        if not str(out.get("reaction_smiles") or "").strip() and _has_concrete_precursors(out):
            p1 = str(out.get("precursor_1") or "").strip()
            p2 = str(out.get("precursor_2") or "").strip()
            out["reaction_smiles"] = f"{p1}.{p2}>>{target_smiles}"
        return out

    inspect = _inspect_target(smiles=target_smiles)
    retrons = _identify_retrons(smiles=target_smiles)
    disconn = _generate_disconnections(smiles=target_smiles, top_k=top_k)

    disconnections = [_with_reaction_smiles(row) for row in list((disconn or {}).get("disconnections", []) or [])]

    # ── Fallback cascade: HTE templates + product similarity ──────────
    hte_template_hits: Dict[str, Any] = {}
    product_similarity_hits: Dict[str, Any] = {}
    fallback_used = False

    if not any(_has_concrete_precursors(row) for row in disconnections):
        # Standard retrons found nothing concrete — try HTE templates.
        hte_template_hits = _apply_hte_templates(target_smiles=target_smiles, top_k=top_k)
        template_disconnections = list(
            (hte_template_hits or {}).get("disconnections", []) or []
        )
        if template_disconnections:
            promoted_rows: List[Dict[str, Any]] = []
            # Promote template hits into the disconnections list
            for rank, td in enumerate(template_disconnections, 1):
                promoted_rows.append({
                    "rank": rank,
                    "reaction_name": td.get("template_name", ""),
                    "retron_name": td.get("template_name", ""),
                    "precursor_1": td.get("precursor_1", ""),
                    "precursor_2": td.get("precursor_2", ""),
                    "reaction_smiles": td.get("reaction_smiles", ""),
                    "description": td.get("description", ""),
                    "difficulty": td.get("difficulty", 0.5),
                    "overall_score": 1.0 - td.get("difficulty", 0.5),
                    "notes": td.get("notes", ""),
                    "source": "hte_template",
                    "hte_conditions": td.get("hte_conditions"),
                })
            disconnections = [_with_reaction_smiles(row) for row in promoted_rows] + disconnections
            fallback_used = True

        if not any(_has_concrete_precursors(row) for row in disconnections):
            # HTE templates also found nothing — try product similarity
            product_similarity_hits = _search_by_product_similarity(
                target_smiles=target_smiles, top_k=top_k,
            )
            sim_precedents = list(
                (product_similarity_hits or {}).get("precedents", []) or []
            )
            if sim_precedents:
                promoted_rows = []
                for rank, sp in enumerate(sim_precedents, 1):
                    promoted_rows.append({
                        "rank": rank,
                        "reaction_name": sp.get("rxn_type", ""),
                        "retron_name": "",
                        "precursor_1": sp.get("precursor_1", ""),
                        "precursor_2": sp.get("precursor_2", ""),
                        "reaction_smiles": sp.get("reaction_smiles", ""),
                        "description": f"Data-driven: {sp.get('product_similarity', 0):.0%} product similarity",
                        "difficulty": round(1.0 - (sp.get("product_similarity") or 0.5), 2),
                        "overall_score": sp.get("product_similarity", 0.5),
                        "notes": f"From HTE {sp.get('rxn_type', 'unknown')} family, yield {sp.get('yield', '?')}%",
                        "source": "product_similarity",
                    })
                disconnections = [_with_reaction_smiles(row) for row in promoted_rows] + disconnections
                fallback_used = True
    else:
        # Retrons succeeded — still run HTE templates in the background
        # for supplementary coverage (cheap call, cached templates)
        hte_template_hits = _apply_hte_templates(target_smiles=target_smiles, top_k=top_k)

    # Keep the facade payload aligned with the locally enriched disconnection rows.
    if disconn.get("success"):
        disconn = dict(disconn)
        disconn["disconnections"] = disconnections
        disconn["total_disconnections"] = len(disconnections)

    # Rebuild disconn dict if fallback injected results
    if fallback_used:
        disconn = {
            "success": True,
            "smiles": target_smiles,
            "total_disconnections": len(disconnections),
            "disconnections": disconnections,
            "fallback_used": True,
        }

    if not disconn.get("success") and not disconnections:
        return _success({
            "target_smiles": target_smiles,
            "inspection": inspect,
            "retrons": retrons,
            "disconnections": disconn,
            "hte_template_hits": hte_template_hits,
            "product_similarity_hits": product_similarity_hits,
        })

    top = next((row for row in disconnections if _has_concrete_precursors(row)), disconnections[0] if disconnections else {})
    top_reaction_name = str(top.get("reaction_name") or "")

    validations: List[Dict[str, Any]] = []
    reaction_smiles_evaluations: List[Dict[str, Any]] = []
    for entry in disconnections[:validate_top_n]:
        p1 = str(entry.get("precursor_1") or "")
        p2 = str(entry.get("precursor_2") or "")
        if not (p1 and p2):
            continue
        validations.append(_check_retro_consistency(
            product_smiles=target_smiles,
            precursor_1=p1,
            precursor_2=p2,
            reaction_name=str(entry.get("reaction_name") or ""),
        ))
        rxn_smiles = str(entry.get("reaction_smiles") or "").strip()
        if rxn_smiles:
            reaction_smiles_evaluations.append(_evaluate_synthesis_proposal(
                mode="reaction",
                reaction_smiles=rxn_smiles,
                reaction_name=str(entry.get("reaction_name") or ""),
            ))

    top_reaction_smiles_evaluation: Dict[str, Any] = {}
    top_rxn_smiles = str(top.get("reaction_smiles") or "").strip()
    if top_rxn_smiles:
        top_reaction_smiles_evaluation = _evaluate_synthesis_proposal(
            mode="reaction",
            reaction_smiles=top_rxn_smiles,
            reaction_name=top_reaction_name,
        )

    conditions: Dict[str, Any] = {}
    if include_conditions and top:
        rxn_smiles = str(top.get("reaction_smiles") or "").strip()
        if not rxn_smiles:
            p1 = str(top.get("precursor_1") or "")
            p2 = str(top.get("precursor_2") or "")
            if p1 and p2:
                rxn_smiles = f"{p1}.{p2}>>{target_smiles}"
        if rxn_smiles:
            conditions = _recommend_conditions_with_strategy(
                reaction_smiles=rxn_smiles,
                top_k=5,
                condition_strategy=condition_strategy,
                condition_source_mode=condition_source_mode,
                condition_selection_mode=condition_selection_mode,
            )

    precedent: Dict[str, Any] = {}
    if include_precedent:
        precedent = _find_retro_precedent(smiles=target_smiles, reaction_name=top_reaction_name)

    result: Dict[str, Any] = {
        "target_smiles": target_smiles,
        "inspection": inspect,
        "retrons": retrons,
        "disconnections": disconn,
        "top_disconnection": top,
        "validation": validations,
        "reaction_smiles_candidates": [
            str(row.get("reaction_smiles") or "").strip()
            for row in disconnections
            if str(row.get("reaction_smiles") or "").strip()
        ],
        "reaction_smiles_evaluations": reaction_smiles_evaluations,
        "top_reaction_smiles_evaluation": top_reaction_smiles_evaluation,
        "conditions_for_top_disconnection": conditions,
        "precedent": precedent,
        "hte_template_hits": hte_template_hits,
        "product_similarity_hits": product_similarity_hits,
        "summary": {
            "total_retrons": int((retrons or {}).get("total_matched", 0) or 0),
            "total_disconnections": len(disconnections),
            "top_reaction_name": top_reaction_name,
            "top_reaction_smiles": top_rxn_smiles,
            "fallback_used": fallback_used,
        },
    }
    return _success(result)


def _forward_synthesis_step(
    smiles_a: str,
    smiles_b: str = "",
    reaction_name: str = "",
    top_k_products: int = 5,
    include_precedent: bool = False,
    include_conditions: bool = True,
    condition_strategy: str = "auto",
    condition_source_mode: str = "",
    condition_selection_mode: str = "best",
) -> Dict[str, Any]:
    """Single-call forward synthesis pass: inspect -> identify -> generate -> rank (+ optional precedent/conditions)."""
    from .forward_synthesis import (
        _find_forward_precedent,
        _generate_products,
        _identify_reactions,
        _inspect_reactants,
        _rank_products,
        _search_reactant_precedent,
    )

    top_k_products = max(1, min(int(top_k_products or 5), 10))

    inspect = _inspect_reactants(smiles_a=smiles_a, smiles_b=smiles_b)
    reactions = _identify_reactions(smiles_a=smiles_a, smiles_b=smiles_b)
    products = _generate_products(
        smiles_a=smiles_a,
        smiles_b=smiles_b,
        reaction_name=reaction_name,
        top_k=top_k_products,
    )
    ranked: Dict[str, Any] = {}
    if products.get("success"):
        ranked = _rank_products(
            products=products.get("products"),
            smiles_a=smiles_a,
            smiles_b=smiles_b,
            top_k=min(3, top_k_products),
        )

    top_product: Dict[str, Any] = {}
    product_list = list((products or {}).get("products", []) or [])
    if product_list:
        top_product = product_list[0]

    product_smiles = str(top_product.get("product_smiles") or "")
    inferred_reaction = (
        str(top_product.get("taxonomy_id") or "")
        or str(top_product.get("template_name") or "")
        or str(reaction_name or "")
    )

    conditions: Dict[str, Any] = {}
    if include_conditions:
        reactants_str = f"{smiles_a}.{smiles_b}" if smiles_b else smiles_a
        if product_smiles:
            rxn_smiles_raw = f"{reactants_str}>>{product_smiles}"
            rxn_smiles, rxn_err = _validate_reaction_smiles(rxn_smiles_raw, require_product=True)
            if rxn_err:
                conditions = _error(f"Invalid top-product reaction SMILES: {rxn_err}")
            else:
                conditions = _recommend_conditions_with_strategy(
                    reaction_smiles=rxn_smiles,
                    top_k=5,
                    condition_strategy=condition_strategy,
                    condition_source_mode=condition_source_mode,
                    condition_selection_mode=condition_selection_mode,
                )
        else:
            conditions = _success({
                "recommendations": [],
                "note": "Skipped condition recommendation: no valid top product was generated.",
            })

    precedent: Dict[str, Any] = {}
    hte_precedent: Dict[str, Any] = {}
    if include_precedent:
        precedent = _find_forward_precedent(
            smiles_a=smiles_a,
            smiles_b=smiles_b,
            reaction_name=reaction_name or inferred_reaction,
        )
        hte_precedent = _search_reactant_precedent(
            smiles_a=smiles_a,
            smiles_b=smiles_b,
            reaction_name=reaction_name or inferred_reaction,
            product_smiles=product_smiles,
            top_k=5,
        )

    return _success({
        "smiles_a": smiles_a,
        "smiles_b": smiles_b,
        "reactant_analysis": inspect,
        "compatible_reactions": reactions,
        "products": products,
        "ranked_products": ranked,
        "top_product": top_product,
        "conditions_for_top_product": conditions,
        "precedent": precedent,
        "hte_precedent": hte_precedent,
        "summary": {
            "total_reactions": int((reactions or {}).get("total_matched", 0) or 0),
            "total_products": int((products or {}).get("total_products", 0) or 0),
            "top_product_smiles": product_smiles,
            "inferred_reaction": inferred_reaction,
        },
    })


def _validate_synthesis_proposal(
    mode: str = "auto",
    reaction_smiles: str = "",
    product_smiles: str = "",
    precursor_1: str = "",
    precursor_2: str = "",
    reaction_name: str = "",
) -> Dict[str, Any]:
    """Backward-compatible validator alias.

    Prefer `evaluate_synthesis_proposal` for richer scoring and consistency checks.
    """
    return _evaluate_synthesis_proposal(
        mode=mode,
        reaction_smiles=reaction_smiles,
        product_smiles=product_smiles,
        precursor_1=precursor_1,
        precursor_2=precursor_2,
        reaction_name=reaction_name,
        include_consistency_checks=True,
    )


def _grade_from_score(score: float) -> str:
    if score >= 85:
        return "A"
    if score >= 70:
        return "B"
    if score >= 55:
        return "C"
    if score >= 40:
        return "D"
    return "F"


def _evaluate_synthesis_proposal(
    mode: str = "auto",
    reaction_smiles: str = "",
    product_smiles: str = "",
    precursor_1: str = "",
    precursor_2: str = "",
    reaction_name: str = "",
    include_consistency_checks: bool = True,
) -> Dict[str, Any]:
    """Dedicated synthesis evaluator with multi-dimensional scoring.

    Dimensions:
      1) Structural validity and chemistry sanity (reaction_eval verdict/checks)
      2) Retrosynthesis complexity sanity (retro mode)
      3) Taxonomy consistency (reaction mode, optional)
    """
    from .reaction_eval import _check_retro_consistency, _evaluate_reaction

    mode_norm = (mode or "auto").strip().lower()
    if mode_norm not in {"auto", "reaction", "retro"}:
        return _error("mode must be one of: auto, reaction, retro")

    if mode_norm == "auto":
        if product_smiles and precursor_1 and precursor_2:
            mode_norm = "retro"
        elif reaction_smiles:
            mode_norm = "reaction"
        else:
            return _error("Provide reaction_smiles or (product_smiles + precursor_1 + precursor_2)")

    if mode_norm == "retro":
        result = _check_retro_consistency(
            product_smiles=product_smiles,
            precursor_1=precursor_1,
            precursor_2=precursor_2,
            reaction_name=reaction_name,
        )
        tool_used = "check_retro_consistency"
    else:
        result = _evaluate_reaction(reaction_smiles=reaction_smiles, reaction_type=reaction_name)
        tool_used = "evaluate_reaction"

    if not isinstance(result, dict):
        return _error("Unexpected evaluation result")
    if not result.get("success"):
        return result

    out = dict(result)
    out["mode_used"] = mode_norm
    out["tool_used"] = tool_used

    # Base score from evaluator verdict
    verdict = str(result.get("verdict") or "").upper()
    verdict_base = {
        "PASS": 90.0,
        "PASS_WITH_WARNINGS": 72.0,
        "FAIL": 35.0,
    }.get(verdict, 50.0)

    critical = len(list(result.get("critical_failures", []) or []))
    warns = len(list(result.get("warnings", []) or []))
    score = verdict_base - (critical * 8.0) - (warns * 2.0)

    consistency: Dict[str, Any] = {
        "checked": False,
        "matches_expected_reaction": None,
        "expected_reaction": reaction_name or "",
        "detected_reaction": "",
        "detected_reaction_id": "",
        "confidence": None,
    }
    if mode_norm == "reaction" and include_consistency_checks:
        consistency["checked"] = True
        try:
            from .chemistry import _detect_reaction_type
            from chemtools.taxonomy.reaction_catalog import resolve_reaction_type

            detected = _detect_reaction_type(reaction_smiles=reaction_smiles)
            if isinstance(detected, dict) and detected.get("success"):
                detected_id = str(detected.get("reaction_type_id") or detected.get("reaction_type") or "")
                consistency["detected_reaction"] = str(detected.get("reaction_type") or "")
                consistency["detected_reaction_id"] = detected_id
                consistency["confidence"] = detected.get("confidence")
                expected_norm = resolve_reaction_type(str(reaction_name)) if reaction_name else ""
                if expected_norm:
                    matches = (expected_norm == detected_id)
                    consistency["matches_expected_reaction"] = bool(matches)
                    if not matches:
                        score -= 15.0
                else:
                    consistency["matches_expected_reaction"] = None
            else:
                score -= 8.0
        except Exception as exc:
            consistency["error"] = str(exc)
            score -= 6.0

    if mode_norm == "retro":
        cx = result.get("complexity_check") if isinstance(result, dict) else {}
        simplified = bool((cx or {}).get("simplification_achieved", True))
        if not simplified:
            score -= 20.0

    score = max(0.0, min(100.0, round(score, 1)))
    out["evaluation_score"] = score
    out["evaluation_grade"] = _grade_from_score(score)
    out["evaluation_summary"] = {
        "verdict": verdict,
        "critical_failures": critical,
        "warnings": warns,
        "consistency": consistency,
    }
    return out


resolve_chemical_tool = ToolPlugin(
    name="resolve_chemical",
    category="composite",
    description=(
        "Unified identifier/SMILES resolver. Use mode='to_smiles' for names/CAS/InChI -> SMILES, "
        "mode='from_smiles' for names/CAS metadata from SMILES, or mode='auto'."
    ),
    prerequisites=[],
    fn=_resolve_chemical,
    provides=["resolved_smiles", "smiles", "canonical_smiles"],
)

reagent_assistant_tool = ToolPlugin(
    name="reagent_assistant",
    category="composite",
    description=(
        "Unified reagent helper. mode='lookup' searches by name/abbreviation/CAS; "
        "mode='list' lists reagents by role (base, ligand, solvent, catalyst) with optional family filter."
    ),
    prerequisites=[],
    fn=_reagent_assistant,
)

search_knowledge_tool = ToolPlugin(
    name="search_knowledge",
    category="composite",
    description=(
        "Search notes and literature together in one call. Optionally include a targeted reaction note "
        "by reaction_note_id (e.g. 'suzuki_miyaura')."
    ),
    prerequisites=[],
    fn=_search_knowledge,
)

read_knowledge_tool = ToolPlugin(
    name="read_knowledge",
    category="composite",
    description=(
        "Read a note, saved literature source, or fetch a URL directly (source='auto'|'notes'|'literature'|'url'). "
        "Use this instead of separate read_notes/read_literature_source/fetch_webpage calls for most cases."
    ),
    prerequisites=[],
    fn=_read_knowledge,
)

analyze_reaction_tool = ToolPlugin(
    name="analyze_reaction",
    category="composite",
    description=(
        "High-level reaction analysis in one call: normalize reaction SMILES, detect taxonomy reaction type, "
        "inspect functional groups, analyze bond changes, and optionally recommend HTE-backed conditions. "
        "For full cross-source condition analysis, set condition_strategy='full' to run literature, motif, "
        "similarity, and rules modes separately, then merge into a consensus-ranked recommendation list. "
        "Use condition_source_mode to force a single source ('literature'|'motif'|'similarity'|'rules')."
    ),
    prerequisites=[],
    fn=_analyze_reaction,
    provides=["reaction_type", "reaction_family", "conditions"],
)

recommend_reaction_conditions_tool = ToolPlugin(
    name="recommend_reaction_conditions",
    category="composite",
    description=(
        "Explicit condition-recommendation facade for reaction SMILES. "
        "Returns ranked catalyst/ligand/base/solvent/temperature options from HTE-backed evidence. "
        "Set condition_strategy='full' to compare literature/motif/similarity/rules and merge consensus."
    ),
    prerequisites=[],
    fn=_recommend_reaction_conditions,
    provides=["conditions", "recommendations"],
)

retrosynthesis_step_tool = ToolPlugin(
    name="retrosynthesis_step",
    category="composite",
    description=(
        "Single-step retrosynthesis facade with automatic fallback cascade: "
        "inspect target → identify retrons → generate disconnections → "
        "(fallback: HTE templates → product similarity when retrons find nothing) → "
        "validate top proposals → suggest conditions and precedent. "
        "Always returns hte_template_hits alongside standard retrons for supplementary coverage. "
        "For robust condition recommendations on the top disconnection, use condition_strategy='full' "
        "to compare literature/motif/similarity/rules and return a merged consensus ranking."
    ),
    prerequisites=[],
    fn=_retrosynthesis_step,
)

forward_synthesis_step_tool = ToolPlugin(
    name="forward_synthesis_step",
    category="composite",
    description=(
        "Single-step forward synthesis facade: inspect reactants, identify compatible reactions, generate and rank products, "
        "and optionally suggest conditions and precedent. "
        "Supports condition_strategy='full' for cross-source condition analysis (literature, motif, "
        "similarity, rules) merged into a consensus-ranked output."
    ),
    prerequisites=[],
    fn=_forward_synthesis_step,
)

validate_synthesis_proposal_tool = ToolPlugin(
    name="validate_synthesis_proposal",
    category="composite",
    description=(
        "Unified reaction/retro validator. Validates reaction SMILES (mode='reaction') or a retrosynthetic "
        "disconnection (mode='retro') using RDKit-based sanity checks and complexity checks."
    ),
    prerequisites=[],
    fn=_validate_synthesis_proposal,
)

evaluate_synthesis_proposal_tool = ToolPlugin(
    name="evaluate_synthesis_proposal",
    category="composite",
    description=(
        "Dedicated synthesis evaluator. Beyond SMILES validity, it combines RDKit reaction sanity checks, "
        "retrosynthesis complexity sanity, optional taxonomy-consistency checks, and returns an overall "
        "evaluation_score (0-100) plus evaluation_grade."
    ),
    prerequisites=[],
    fn=_evaluate_synthesis_proposal,
)


def _project_analyze_reaction(result: dict) -> Dict[str, Any]:
    """Project nested analyze_reaction outputs into top-level structured fields."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    out: Dict[str, Any] = {}
    det = result.get("reaction_type") or {}
    if isinstance(det, dict) and det.get("success"):
        out["reaction_type"] = det.get("reaction_type_id") or det.get("reaction_type")
        out["reaction_family"] = det.get("family_label")
        if det.get("reaction_type_metadata"):
            out["reaction_type_metadata"] = det.get("reaction_type_metadata")
    bchg = result.get("bond_changes") or {}
    if isinstance(bchg, dict) and bchg.get("success"):
        out["bonds_formed"] = bchg.get("bonds_formed", [])
        out["bonds_broken"] = bchg.get("bonds_broken", [])
        out["key_bond_type"] = bchg.get("key_bond_type")
    cond = result.get("conditions") or {}
    if isinstance(cond, dict) and cond.get("success"):
        out["conditions"] = cond.get("recommendations", [])
    return {k: v for k, v in out.items() if v is not None}


def _project_recommend_reaction_conditions(result: dict) -> Dict[str, Any]:
    """Project recommendation output into structured conditions."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}
    return {
        "conditions": result.get("recommendations", []),
    }


def _project_retrosynthesis_step(result: dict) -> Dict[str, Any]:
    """Project top retrosynthesis route details into structured output."""
    if not isinstance(result, dict) or not result.get("success"):
        return {}

    out: Dict[str, Any] = {}
    top = result.get("top_disconnection") or {}
    if isinstance(top, dict):
        if top.get("reaction_smiles"):
            out["reaction_smiles"] = top.get("reaction_smiles")
        if top.get("reaction_name"):
            out["reaction_type"] = top.get("reaction_name")
        if top.get("precursor_1") and top.get("precursor_2"):
            out["top_precursors"] = [
                top.get("precursor_1"),
                top.get("precursor_2"),
            ]

    top_eval = result.get("top_reaction_smiles_evaluation") or {}
    if isinstance(top_eval, dict) and top_eval.get("success"):
        out["reaction_smiles_evaluation"] = {
            "verdict": top_eval.get("verdict"),
            "evaluation_score": top_eval.get("evaluation_score"),
            "evaluation_grade": top_eval.get("evaluation_grade"),
            "summary": top_eval.get("evaluation_summary"),
        }

    cond = result.get("conditions_for_top_disconnection") or {}
    if isinstance(cond, dict) and cond.get("success"):
        out["conditions"] = cond.get("recommendations", [])

    return {k: v for k, v in out.items() if v is not None}


analyze_reaction_tool.structured_projection = _project_analyze_reaction
recommend_reaction_conditions_tool.structured_projection = _project_recommend_reaction_conditions
retrosynthesis_step_tool.structured_projection = _project_retrosynthesis_step


COMPOSITE_TOOLS = [
    resolve_chemical_tool,
    reagent_assistant_tool,
    search_knowledge_tool,
    read_knowledge_tool,
    analyze_reaction_tool,
    recommend_reaction_conditions_tool,
    retrosynthesis_step_tool,
    forward_synthesis_step_tool,
    evaluate_synthesis_proposal_tool,
    validate_synthesis_proposal_tool,
]
