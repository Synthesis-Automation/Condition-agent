"""
Detection validation using reacted motifs patterns.

This module provides second-pass validation to correct common misclassifications
using taxonomy constraints from either:
1) reacted/formed motif sets, or
2) CRK reaction keys (preferred path).

NOW TAXONOMY-DRIVEN: Validation patterns are loaded from taxonomy files
(reaction_types.v4.0.json + compound_logic.json) instead of hardcoded.
"""

import re
from typing import Dict, Any, Set, List, Optional, Tuple
from chemtools.taxonomy.reaction_catalog import (
    load_reaction_catalog,
    ReactionTypeDefinition,
)
from .utils import normalize_motif_id


# Cache the taxonomy for performance
_CATALOG_CACHE: Optional[Tuple[Dict[str, ReactionTypeDefinition], Dict[str, str]]] = None


def _get_catalog() -> Tuple[Dict[str, ReactionTypeDefinition], Dict[str, str]]:
    """Load and cache the reaction catalog."""
    global _CATALOG_CACHE
    if _CATALOG_CACHE is None:
        _CATALOG_CACHE = load_reaction_catalog()
    return _CATALOG_CACHE


def _motifs_match_slot(motifs: Set[str], allowed: List[str]) -> bool:
    """Check if any motif matches the allowed patterns for a slot."""
    return any(m in motifs for m in allowed)


def _as_str_set(values: Any) -> Set[str]:
    if not values:
        return set()
    if isinstance(values, (list, tuple, set)):
        return {str(v) for v in values if str(v).strip()}
    return {str(values)} if str(values).strip() else set()


def _constraints_match(
    constraints: Dict[str, Any],
    reacted_set: Set[str],
    formed_set: Set[str],
    *,
    formed_bond_tokens: Optional[Set[str]] = None,
    broken_bond_tokens: Optional[Set[str]] = None,
    reactant_slot_matches: int = 0,
    product_slot_matches: int = 0,
) -> bool:
    """Apply lightweight include/exclude constraints on reacted/formed motifs."""
    include_reacted = _as_str_set(constraints.get("include_reacted"))
    exclude_reacted = _as_str_set(constraints.get("exclude_reacted"))
    include_formed = _as_str_set(constraints.get("include_formed"))
    exclude_formed = _as_str_set(constraints.get("exclude_formed"))
    include_bond_formed = _as_str_set(constraints.get("include_bond_formed"))
    exclude_bond_formed = _as_str_set(constraints.get("exclude_bond_formed"))
    include_bond_broken = _as_str_set(constraints.get("include_bond_broken"))
    exclude_bond_broken = _as_str_set(constraints.get("exclude_bond_broken"))
    min_reactant_slot_matches = int(constraints.get("min_reactant_slot_matches") or 0)
    min_product_slot_matches = int(constraints.get("min_product_slot_matches") or 0)
    formed_bonds = formed_bond_tokens or set()
    broken_bonds = broken_bond_tokens or set()

    if include_reacted and not include_reacted.issubset(reacted_set):
        return False
    if exclude_reacted and (exclude_reacted & reacted_set):
        return False
    if include_formed and not include_formed.issubset(formed_set):
        return False
    if exclude_formed and (exclude_formed & formed_set):
        return False
    # Enforce include-bond constraints only when bond evidence is available.
    if include_bond_formed and formed_bond_tokens and not include_bond_formed.issubset(formed_bonds):
        return False
    if exclude_bond_formed and (exclude_bond_formed & formed_bonds):
        return False
    if include_bond_broken and broken_bond_tokens and not include_bond_broken.issubset(broken_bonds):
        return False
    if exclude_bond_broken and (exclude_bond_broken & broken_bonds):
        return False
    if reactant_slot_matches < max(0, min_reactant_slot_matches):
        return False
    if product_slot_matches < max(0, min_product_slot_matches):
        return False
    return True


_SNAR_ACTIVATING_TOKENS = (
    "HeteroAr",
    "AromN",
    "Pyridine",
    "Pyridyl",
    "Pyrimidine",
    "Pyrimidyl",
    "Pyrazine",
    "Triazine",
    "Quinoline",
    "Isoquinoline",
    "-NO2",
    "-CN",
    "-COR",
    "-CO2",
    "-SO2",
    "-CHO",
    "-CF3",
    "-N2+",
)
_SNAR_DEACTIVATING_TOKENS = (
    "Ar-OR",
    "Ar-NR2",
    "Ar-NHR",
    "Ar-NH2",
    "Ar-SR",
    "Ar-Alkyl",
)


def _snar_activation_supported(
    reacted_set: Set[str],
    formed_set: Set[str],
    spectators_set: Optional[Set[str]] = None,
) -> bool:
    """Return whether motif evidence supports SNAr electronic activation."""
    motif_space: Set[str] = set(reacted_set) | set(formed_set) | set(spectators_set or set())
    if not motif_space:
        return False

    # Heteroaryl context is a strong SNAr activation signal.
    for motif in motif_space:
        text = str(motif)
        if text.startswith("HeteroAr-") or text.startswith("AromN-"):
            return True
        if any(tag in text for tag in ("Pyridine", "Pyridyl", "Pyrimidine", "Pyrazine", "Triazine", "Quinoline")):
            return True

    # Common electron-withdrawing activation tags.
    for motif in motif_space:
        text = str(motif)
        if any(token in text for token in _SNAR_ACTIVATING_TOKENS):
            return True

    # Explicitly down-rank electron-rich aryl systems without activation.
    has_deactivating = any(
        any(token in str(motif) for token in _SNAR_DEACTIVATING_TOKENS)
        for motif in motif_space
    )
    if has_deactivating:
        return False

    # Plain Ar-X without activation defaults to unsupported SNAr assignment.
    return False


def _score_specificity_candidate(
    *,
    reaction_id: str,
    matched_slots: List[str],
    product_slot_matches: int,
    reactant_support: int,
    product_support: int,
    reactant_allowed_total: int,
    product_allowed_total: int,
    catalog_index: int,
) -> Tuple[int, int, int, int, int, int, int, int]:
    """Score candidate specificity for stable ranking."""
    return (
        len(matched_slots),          # prefer broader slot support
        product_slot_matches,        # prefer stronger product evidence
        reactant_support,            # reactant evidence overlap
        product_support,             # product evidence overlap
        -reactant_allowed_total,     # fewer allowed motifs => more specific
        -product_allowed_total,      # fewer allowed motifs => more specific
        -len(reaction_id),           # stable tie-break
        -catalog_index,              # preserve deterministic order as last tie-break
    )


def _split_motif_tokens(value: str) -> List[str]:
    tokens = [tok.strip() for tok in re.split(r"[|,;/]", value) if tok.strip()]
    return [normalize_motif_id(tok) for tok in tokens if tok and tok != "[]"]


def _canonicalize_bond_token(token: str) -> Optional[str]:
    text = str(token or "").strip()
    if "-" not in text:
        return None
    left, right = [t.strip() for t in text.split("-", 1)]
    left_match = re.search(r"[A-Z][a-z]?", left)
    right_match = re.search(r"[A-Z][a-z]?", right)
    if not left_match or not right_match:
        return None
    left_el = left_match.group(0)
    right_el = right_match.group(0)
    return "-".join(sorted((left_el, right_el)))


def _parse_crk_key(
    reaction_key: str,
) -> Tuple[Set[str], Set[str], Set[str], List[str], List[str]]:
    """
    Parse CRK key into reacted/formed/spectator motif sets and bond tokens.

    Returns:
        (reacted_set, formed_set, spectators_set, bond_formed_tokens, bond_broken_tokens)
    """
    if not reaction_key:
        return set(), set(), set(), [], []

    text = str(reaction_key).strip()
    if not text or "->" not in text:
        return set(), set(), set(), [], []

    sections = [section.strip() for section in text.split(" | ") if section.strip()]
    if not sections:
        return set(), set(), set(), [], []

    summary = sections[0]
    if summary.startswith("CRK-v1"):
        summary = summary[len("CRK-v1"):].strip()
    if summary.startswith("|"):
        summary = summary[1:].strip()

    reacted: Set[str] = set()
    formed: Set[str] = set()
    if "->" in summary:
        reactant_part, product_part = summary.split("->", 1)
        reacted = set(_split_motif_tokens(reactant_part.strip()))
        formed = set(_split_motif_tokens(product_part.strip()))

    spectators: Set[str] = set()
    formed_bonds: List[str] = []
    broken_bonds: List[str] = []
    for section in sections[1:]:
        lower = section.lower()
        if lower.startswith("spectators:"):
            payload = section.split(":", 1)[1].strip()
            spectators = set(_split_motif_tokens(payload))
        elif lower.startswith("bond_formed:"):
            payload = section.split(":", 1)[1].strip()
            formed_bonds = [tok.strip() for tok in payload.split(";") if tok.strip()]
        elif lower.startswith("bond_broken:"):
            payload = section.split(":", 1)[1].strip()
            broken_bonds = [tok.strip() for tok in payload.split(";") if tok.strip()]

    return reacted, formed, spectators, formed_bonds, broken_bonds


def _collect_ranked_catalog_candidates(
    reacted_set: Set[str],
    formed_set: Set[str],
    spectators_set: Optional[Set[str]] = None,
    *,
    formed_bond_tokens: Optional[Set[str]] = None,
    broken_bond_tokens: Optional[Set[str]] = None,
) -> List[Dict[str, Any]]:
    """Collect and rank taxonomy candidates using specificity-aware scoring."""
    definitions, _ = _get_catalog()
    candidates: List[Dict[str, Any]] = []

    for catalog_index, (reaction_id, defn) in enumerate(definitions.items()):
        has_reactant_slots = bool(defn.reactants)
        has_product_slots = bool(defn.products)
        if not has_reactant_slots and not has_product_slots:
            continue

        all_slots_match = True
        matched_slots: List[str] = []
        reactant_allowed_union: Set[str] = set()
        reactant_allowed_total = 0
        for slot_name, slot_req in defn.reactants.items():
            if not slot_req.allowed:
                continue
            allowed = set(slot_req.allowed)
            reactant_allowed_union.update(allowed)
            reactant_allowed_total += len(allowed)
            if _motifs_match_slot(reacted_set, slot_req.allowed):
                matched_slots.append(slot_name)
            else:
                all_slots_match = False
                break

        if not all_slots_match:
            continue
        if has_reactant_slots and len(matched_slots) < 1:
            continue

        product_match = True
        product_slot_matches = 0
        product_allowed_union: Set[str] = set()
        product_allowed_total = 0
        if defn.products:
            product_match = False
            for slot_req in defn.products.values():
                if slot_req.allowed and _motifs_match_slot(formed_set, slot_req.allowed):
                    product_match = True
                    product_slot_matches += 1
                if slot_req.allowed:
                    allowed = set(slot_req.allowed)
                    product_allowed_union.update(allowed)
                    product_allowed_total += len(allowed)

        if not product_match:
            continue

        if str(reaction_id) == "SNAr_CN":
            if not _snar_activation_supported(
                reacted_set,
                formed_set,
                spectators_set=spectators_set,
            ):
                continue

        if not _constraints_match(
            defn.constraints or {},
            reacted_set,
            formed_set,
            formed_bond_tokens=formed_bond_tokens,
            broken_bond_tokens=broken_bond_tokens,
            reactant_slot_matches=len(matched_slots),
            product_slot_matches=product_slot_matches,
        ):
            continue

        reactant_support = len(reacted_set & reactant_allowed_union)
        product_support = len(formed_set & product_allowed_union)
        score = _score_specificity_candidate(
            reaction_id=reaction_id,
            matched_slots=matched_slots,
            product_slot_matches=product_slot_matches,
            reactant_support=reactant_support,
            product_support=product_support,
            reactant_allowed_total=reactant_allowed_total,
            product_allowed_total=product_allowed_total,
            catalog_index=catalog_index,
        )

        candidates.append(
            {
                "reaction_id": reaction_id,
                "display_name": defn.name,
                "matched_reactant_slots": list(matched_slots),
                "matched_product_slots": int(product_slot_matches),
                "reactant_support": int(reactant_support),
                "product_support": int(product_support),
                "reactant_allowed_total": int(reactant_allowed_total),
                "product_allowed_total": int(product_allowed_total),
                "score": score,
                "catalog_index": int(catalog_index),
            }
        )

    candidates.sort(key=lambda c: c.get("score"), reverse=True)
    return candidates


def _match_reaction_catalog_specificity(
    reacted_set: Set[str],
    formed_set: Set[str],
    spectators_set: Optional[Set[str]] = None,
    *,
    formed_bond_tokens: Optional[Set[str]] = None,
    broken_bond_tokens: Optional[Set[str]] = None,
) -> Optional[Tuple[str, List[str], str]]:
    """Return best match using specificity-aware ranked candidates."""
    candidates = _collect_ranked_catalog_candidates(
        reacted_set,
        formed_set,
        spectators_set=spectators_set,
        formed_bond_tokens=formed_bond_tokens,
        broken_bond_tokens=broken_bond_tokens,
    )
    if not candidates:
        return None
    top = candidates[0]
    return (
        str(top["reaction_id"]),
        [str(s) for s in top["matched_reactant_slots"]],
        str(top["display_name"]),
    )


def _match_reaction_catalog_legacy(
    reacted_set: Set[str],
    formed_set: Set[str],
    spectators_set: Optional[Set[str]] = None,
    *,
    formed_bond_tokens: Optional[Set[str]] = None,
    broken_bond_tokens: Optional[Set[str]] = None,
) -> Optional[Tuple[str, List[str], str]]:
    """Return first taxonomy match (legacy pre-specificity behavior)."""
    definitions, _ = _get_catalog()
    for reaction_id, defn in definitions.items():
        has_reactant_slots = bool(defn.reactants)
        has_product_slots = bool(defn.products)
        if not has_reactant_slots and not has_product_slots:
            continue

        all_slots_match = True
        matched_slots: List[str] = []
        for slot_name, slot_req in defn.reactants.items():
            if not slot_req.allowed:
                continue
            if _motifs_match_slot(reacted_set, slot_req.allowed):
                matched_slots.append(slot_name)
            else:
                all_slots_match = False
                break

        if not all_slots_match:
            continue
        if has_reactant_slots and len(matched_slots) < 1:
            continue

        product_match = True
        if defn.products:
            product_match = False
            for slot_req in defn.products.values():
                if slot_req.allowed and _motifs_match_slot(formed_set, slot_req.allowed):
                    product_match = True
                    break
        if not product_match:
            continue
        if str(reaction_id) == "SNAr_CN":
            if not _snar_activation_supported(
                reacted_set,
                formed_set,
                spectators_set=spectators_set,
            ):
                continue
        if not _constraints_match(
            defn.constraints or {},
            reacted_set,
            formed_set,
            formed_bond_tokens=formed_bond_tokens,
            broken_bond_tokens=broken_bond_tokens,
            reactant_slot_matches=len(matched_slots),
            product_slot_matches=1 if product_match else 0,
        ):
            continue
        return reaction_id, matched_slots, defn.name
    return None


def _match_reaction_catalog(
    reacted_set: Set[str],
    formed_set: Set[str],
    spectators_set: Optional[Set[str]] = None,
    *,
    use_legacy: bool = False,
    formed_bond_tokens: Optional[Set[str]] = None,
    broken_bond_tokens: Optional[Set[str]] = None,
) -> Optional[Tuple[str, List[str], str]]:
    """Return best taxonomy match as (reaction_id, matched_slots, display_name)."""
    if use_legacy:
        return _match_reaction_catalog_legacy(
            reacted_set,
            formed_set,
            spectators_set=spectators_set,
            formed_bond_tokens=formed_bond_tokens,
            broken_bond_tokens=broken_bond_tokens,
        )
    return _match_reaction_catalog_specificity(
        reacted_set,
        formed_set,
        spectators_set=spectators_set,
        formed_bond_tokens=formed_bond_tokens,
        broken_bond_tokens=broken_bond_tokens,
    )


def _build_match_evidence(
    *,
    reacted_set: Set[str],
    formed_set: Set[str],
    spectators_set: Optional[Set[str]] = None,
    selected_match: Optional[Tuple[str, List[str], str]],
    use_legacy: bool,
    max_candidates: int = 5,
) -> Dict[str, Any]:
    matcher = "taxonomy_legacy_v1" if use_legacy else "taxonomy_specificity_v2"
    candidates = _collect_ranked_catalog_candidates(
        reacted_set,
        formed_set,
        spectators_set=spectators_set,
    )
    top_candidates = [
        {
            "reaction_id": c.get("reaction_id"),
            "display_name": c.get("display_name"),
            "matched_reactant_slots": c.get("matched_reactant_slots"),
            "matched_product_slots": c.get("matched_product_slots"),
            "reactant_support": c.get("reactant_support"),
            "product_support": c.get("product_support"),
            "reactant_allowed_total": c.get("reactant_allowed_total"),
            "product_allowed_total": c.get("product_allowed_total"),
            "score": list(c.get("score") or []),
        }
        for c in candidates[: max(1, int(max_candidates))]
    ]

    selected_payload: Optional[Dict[str, Any]] = None
    if selected_match:
        rid, matched_slots, display_name = selected_match
        selected_payload = {
            "reaction_id": rid,
            "display_name": display_name,
            "matched_reactant_slots": list(matched_slots),
        }
        candidate = next((c for c in candidates if c.get("reaction_id") == rid), None)
        if candidate:
            selected_payload["matched_product_slots"] = candidate.get("matched_product_slots")
            selected_payload["score"] = list(candidate.get("score") or [])

    return {
        "matcher": matcher,
        "reacted_motifs": sorted(reacted_set),
        "formed_motifs": sorted(formed_set),
        "selected": selected_payload,
        "top_candidates": top_candidates,
    }


def _format_validation_response(
    *,
    initial_detection: str,
    initial_confidence: float,
    match: Optional[Tuple[str, List[str], str]],
    method: str,
    reason_prefix: str,
    evidence: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    if match is None:
        payload = {
            "reaction_type": initial_detection,
            "confidence": initial_confidence,
            "validation_method": "slot_based_confirmed",
            "corrected_from": None,
            "reason": "Pattern consistent with slot-based detection",
        }
        initial_text = str(initial_detection or "").strip()
        if initial_text and initial_text != "Unknown" and isinstance(evidence, dict):
            top_candidates = evidence.get("top_candidates")
            reacted_evidence = evidence.get("reacted_motifs") or []
            formed_evidence = evidence.get("formed_motifs") or []
            if (
                isinstance(top_candidates, list)
                and not top_candidates
                and (reacted_evidence or formed_evidence)
            ):
                payload = {
                    "reaction_type": "Unknown",
                    "confidence": 0.0,
                    "validation_method": "crk_pattern_conflict",
                    "corrected_from": initial_detection,
                    "reason": "No taxonomy-consistent CRK match for observed motif/bond evidence",
                }
        if evidence is not None:
            payload["evidence"] = evidence
        return payload

    reaction_id, matched_slots, display_name = match
    if reaction_id != initial_detection:
        payload = {
            "reaction_type": reaction_id,
            "confidence": 0.95,
            "validation_method": method,
            "corrected_from": initial_detection,
            "reason": f"{reason_prefix}: {' + '.join(matched_slots)} -> {display_name}",
        }
        if evidence is not None:
            payload["evidence"] = evidence
        return payload

    payload = {
        "reaction_type": initial_detection,
        "confidence": max(initial_confidence, 0.95),
        "validation_method": "slot_based_confirmed",
        "corrected_from": None,
        "reason": "Pattern consistent with slot-based detection",
    }
    if evidence is not None:
        payload["evidence"] = evidence
    return payload


def validate_detection_with_crk_key(
    initial_detection: str,
    initial_confidence: float,
    reaction_key: str,
    *,
    use_legacy: bool = False,
    include_evidence: bool = True,
    max_candidates: int = 5,
) -> Dict[str, Any]:
    """
    Validate and refine reaction type detection from CRK key.

    This is the preferred streamlined path: CRK_raw -> taxonomy match.
    """
    reacted_set, formed_set, _spectators, _formed_bonds, _broken_bonds = _parse_crk_key(reaction_key)
    formed_bond_tokens = {
        tok
        for tok in (_canonicalize_bond_token(t) for t in (_formed_bonds or []))
        if tok
    }
    broken_bond_tokens = {
        tok
        for tok in (_canonicalize_bond_token(t) for t in (_broken_bonds or []))
        if tok
    }
    match = _match_reaction_catalog(
        reacted_set,
        formed_set,
        spectators_set=_spectators,
        use_legacy=use_legacy,
        formed_bond_tokens=formed_bond_tokens,
        broken_bond_tokens=broken_bond_tokens,
    )
    evidence = (
        _build_match_evidence(
            reacted_set=reacted_set,
            formed_set=formed_set,
            spectators_set=_spectators,
            selected_match=match,
            use_legacy=use_legacy,
            max_candidates=max_candidates,
        )
        if include_evidence
        else None
    )
    return _format_validation_response(
        initial_detection=initial_detection,
        initial_confidence=initial_confidence,
        match=match,
        method="crk_pattern",
        reason_prefix="Taxonomy pattern (CRK)",
        evidence=evidence,
    )


def validate_detection_with_reacted_motifs(
    initial_detection: str,
    initial_confidence: float,
    reacted_motifs: List[str],
    formed_motifs: List[str],
    spectator_motifs: Optional[List[str]] = None,
    *,
    use_legacy: bool = False,
    include_evidence: bool = True,
    max_candidates: int = 5,
) -> Dict[str, Any]:
    """
    Validate and refine reaction type detection using reacted motifs.
    
    This adds a logical second-pass that corrects common misclassifications
    by using actual motif consumption patterns.
    
    NOW TAXONOMY-DRIVEN: Loads patterns from reaction_types.v4.0.json + compound_logic.json
    instead of hardcoding motif lists.
    
    Args:
        initial_detection: Reaction type from slot-based detection
        initial_confidence: Confidence from slot-based detection
        reacted_motifs: Motifs consumed in reaction (from aggregates)
        formed_motifs: Motifs formed in products (from aggregates)
        spectator_motifs: Motifs present but unchanged (optional)
        
    Returns:
        Dict with validated reaction_type, confidence, and validation metadata
    """
    reacted_set = {normalize_motif_id(str(m)) for m in (reacted_motifs or []) if m}
    formed_set = {normalize_motif_id(str(m)) for m in (formed_motifs or []) if m}
    spectator_set = {normalize_motif_id(str(m)) for m in (spectator_motifs or []) if m}
    match = _match_reaction_catalog(
        reacted_set,
        formed_set,
        spectators_set=spectator_set,
        use_legacy=use_legacy,
    )
    evidence = (
        _build_match_evidence(
            reacted_set=reacted_set,
            formed_set=formed_set,
            spectators_set=spectator_set,
            selected_match=match,
            use_legacy=use_legacy,
            max_candidates=max_candidates,
        )
        if include_evidence
        else None
    )
    return _format_validation_response(
        initial_detection=initial_detection,
        initial_confidence=initial_confidence,
        match=match,
        method="reacted_motifs_pattern",
        reason_prefix="Taxonomy pattern",
        evidence=evidence,
    )
