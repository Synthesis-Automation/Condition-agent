"""
Deterministic reaction event extraction and reaction-key quality assessment.

This module provides a general-purpose, taxonomy-aligned evidence layer that
summarizes what transformation primitives occurred (e.g., C-N formation,
leaving-group displacement, ring closure) and flags anomalies suggesting
multi-step records or missing activation metadata.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.util import rdkit_helpers


_HALOGENS = {"F", "Cl", "Br", "I"}


def _split_reaction_sides(reaction_smiles: str) -> Tuple[List[str], List[str]]:
    if not reaction_smiles:
        return [], []
    text = str(reaction_smiles).strip()
    if ">>" in text:
        left, right = text.split(">>", 1)
        return [x for x in left.split(".") if x], [x for x in right.split(".") if x]
    parts = text.split(">")
    if len(parts) == 3:
        left = parts[0]
        right = parts[2]
        return [x for x in left.split(".") if x], [x for x in right.split(".") if x]
    return [x for x in text.split(".") if x], []


def _extract_bond_section(bond_key: Optional[str], section: str) -> List[str]:
    if not bond_key:
        return []
    marker = f"{section}: "
    parts = [p.strip() for p in str(bond_key).split(" | ") if p.strip()]
    for part in parts:
        if part.startswith(marker):
            payload = part[len(marker):].strip()
            if not payload:
                return []
            return [tok.strip() for tok in payload.split(";") if tok.strip()]
    return []


def _parse_element(label: str) -> str:
    text = str(label).strip()
    if not text:
        return ""
    if text.startswith("#"):
        return ""
    if "(" in text:
        text = text.split("(", 1)[0]
    if "-" in text and len(text) > 2:
        # tokens like "C(ar)-N" are split earlier; this is a safety guard.
        text = text.split("-", 1)[0]
    return text


def _parse_bond_token(token: str) -> Optional[Tuple[str, str]]:
    if "-" not in str(token):
        return None
    left, right = [t.strip() for t in str(token).split("-", 1)]
    left_el = _parse_element(left)
    right_el = _parse_element(right)
    if not left_el or not right_el:
        return None
    pair = tuple(sorted((left_el, right_el)))
    return pair  # type: ignore[return-value]


def _bond_pairs(tokens: Iterable[str]) -> Set[Tuple[str, str]]:
    pairs: Set[Tuple[str, str]] = set()
    for token in tokens or []:
        parsed = _parse_bond_token(token)
        if parsed is not None:
            pairs.add(parsed)
    return pairs


def _ring_count(smiles_list: Iterable[str]) -> int:
    if not rdkit_helpers.rdkit_available():
        return 0
    count = 0
    for smiles in smiles_list or []:
        mol = rdkit_helpers.parse_smiles(smiles)
        if mol is None:
            continue
        try:
            count += int(mol.GetRingInfo().NumRings())
        except Exception:
            continue
    return count


def _has_carboxylic_acid_like(motif_ids: Set[str]) -> bool:
    for mid in motif_ids:
        m = str(mid)
        if "CO2H" in m or "COOH" in m:
            return True
    return False


def _has_explicit_acyl_activation_like(motif_ids: Set[str]) -> bool:
    markers = ("COCl", "COBr", "COF", "Anhydride", "MixedAnhydride", "Acyl-Imidazole")
    for mid in motif_ids:
        m = str(mid)
        if any(marker in m for marker in markers):
            return True
    return False


def _event(kind: str, confidence: float, details: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    payload: Dict[str, Any] = {"kind": kind, "confidence": max(0.0, min(1.0, float(confidence)))}
    if details:
        payload["details"] = details
    return payload


def _quality_bucket(score: float) -> str:
    if score >= 0.75:
        return "high"
    if score >= 0.45:
        return "medium"
    return "low"


def summarize_reaction_events(
    *,
    reaction_smiles: str,
    bond_key: Optional[str],
    fallback_bond_key: Optional[str],
    reacted_motifs: Iterable[str],
    formed_motifs: Iterable[str],
    mapping_warning: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Build deterministic reaction-event summary and key quality diagnostics.
    """
    active_bond_key = bond_key or fallback_bond_key
    formed_tokens = _extract_bond_section(active_bond_key, section="form")
    broken_tokens = _extract_bond_section(active_bond_key, section="break")
    formed_pairs = _bond_pairs(formed_tokens)
    broken_pairs = _bond_pairs(broken_tokens)

    reactant_smiles, product_smiles = _split_reaction_sides(reaction_smiles)
    reactant_ring_count = _ring_count(reactant_smiles)
    product_ring_count = _ring_count(product_smiles)
    ring_delta = product_ring_count - reactant_ring_count

    reacted_set = {str(m) for m in reacted_motifs or [] if str(m).strip()}
    formed_set = {str(m) for m in formed_motifs or [] if str(m).strip()}

    events: List[Dict[str, Any]] = []
    anomalies: List[str] = []

    if ("C", "N") in formed_pairs:
        events.append(_event("c_n_bond_formation", 0.9))
    if ("C", "O") in formed_pairs:
        events.append(_event("c_o_bond_formation", 0.9))
    if ("C", "S") in formed_pairs:
        events.append(_event("c_s_bond_formation", 0.9))
    if ("C", "C") in formed_pairs:
        events.append(_event("c_c_bond_formation", 0.85))

    for hal in _HALOGENS:
        if tuple(sorted(("C", hal))) in broken_pairs:
            nucleophile = None
            for cand in ("N", "O", "S", "C"):
                if tuple(sorted(("C", cand))) in formed_pairs:
                    nucleophile = cand
                    break
            events.append(
                _event(
                    "leaving_group_displacement",
                    0.92,
                    {"leaving_group": hal, "nucleophile_element": nucleophile},
                )
            )

    if ring_delta > 0:
        events.append(
            _event(
                "ring_closure_or_annulation",
                0.88,
                {
                    "reactant_ring_count": reactant_ring_count,
                    "product_ring_count": product_ring_count,
                    "ring_delta": ring_delta,
                },
            )
        )

    if len(reactant_smiles) == 1 and formed_pairs:
        events.append(_event("intramolecular_likely", 0.8))
    elif len(reactant_smiles) > 1 and formed_pairs:
        events.append(_event("intermolecular_or_multi_component", 0.7))

    if ("C", "N") in formed_pairs and _has_carboxylic_acid_like(reacted_set):
        events.append(_event("amidation_like", 0.78))
        if not _has_explicit_acyl_activation_like(reacted_set):
            anomalies.append("amidation_without_explicit_activation_marker")

    primary_kinds = {
        "c_n_bond_formation",
        "c_o_bond_formation",
        "c_s_bond_formation",
        "c_c_bond_formation",
        "leaving_group_displacement",
        "ring_closure_or_annulation",
        "amidation_like",
    }
    primary_count = sum(1 for ev in events if ev.get("kind") in primary_kinds)
    if primary_count >= 3:
        anomalies.append("multi_transform_record_possible")

    if isinstance(mapping_warning, dict) and mapping_warning:
        anomalies.append("mapping_unreliable_fallback_used")

    quality_score = 1.0
    quality_reasons: List[str] = []
    if not formed_tokens and not formed_set:
        quality_score -= 0.45
        quality_reasons.append("missing_formed_bond_and_product_motif_evidence")
    if not active_bond_key:
        quality_score -= 0.25
        quality_reasons.append("missing_bond_key")
    if anomalies:
        quality_score -= min(0.45, 0.15 * len(anomalies))
        quality_reasons.extend(anomalies)
    quality_score = max(0.0, min(1.0, quality_score))

    return {
        "events": events,
        "anomalies": anomalies,
        "bond_pairs": {
            "formed": sorted(list(formed_pairs)),
            "broken": sorted(list(broken_pairs)),
        },
        "ring_change": {
            "reactants": reactant_ring_count,
            "products": product_ring_count,
            "delta": ring_delta,
        },
        "reaction_key_quality": {
            "score_0_1": round(quality_score, 3),
            "level": _quality_bucket(quality_score),
            "reasons": quality_reasons,
        },
    }

