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
_PATTERN_SMARTS = {
    "ester": "[CX3](=O)[OX2][#6]",
    "carboxy_acid": "[CX3](=O)[OX2H]",
    "carboxylate": "[CX3](=O)[O-]",
    "phenol": "[c][OX2H]",
    "benzyl_halide": "[c][CH2][Cl,Br,I]",
    "benzyl_aryl_ether": "[c][CH2][O][c]",
}
_EVENT_SIGNATURE_PRIORITY = [
    "benzyl_o_alkylation_like",
    "ester_hydrolysis_like",
    "amidation_like",
    "ring_closure_or_annulation",
    "leaving_group_displacement",
    "c_n_bond_formation",
    "c_o_bond_formation",
    "c_s_bond_formation",
    "c_c_bond_formation",
]
_EVENT_SIGNATURE_CODE = {
    "benzyl_o_alkylation_like": "BzOAlk",
    "ester_hydrolysis_like": "EsterHyd",
    "amidation_like": "Amid",
    "ring_closure_or_annulation": "Ann",
    "leaving_group_displacement": "LGDisp",
    "c_n_bond_formation": "C-N",
    "c_o_bond_formation": "C-O",
    "c_s_bond_formation": "C-S",
    "c_c_bond_formation": "C-C",
}
_REDOX_HETERO = {"N", "O", "S", "P", "F", "Cl", "Br", "I"}


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


def _compute_redox_assessment(
    formed_pairs: Set[Tuple[str, str]],
    broken_pairs: Set[Tuple[str, str]],
) -> Dict[str, Any]:
    """
    Heuristic redox assessment from bond-pair deltas.

    Score convention:
    - positive: net oxidation signal on substrate carbons
    - negative: net reduction signal on substrate carbons
    """

    signals: List[str] = []
    score = 0.0

    formed_het = sum(
        1 for a, b in formed_pairs if "C" in (a, b) and ({a, b} - {"C"}) and next(iter({a, b} - {"C"})) in _REDOX_HETERO
    )
    broken_het = sum(
        1 for a, b in broken_pairs if "C" in (a, b) and ({a, b} - {"C"}) and next(iter({a, b} - {"C"})) in _REDOX_HETERO
    )
    formed_ch = int(("C", "H") in formed_pairs)
    broken_ch = int(("C", "H") in broken_pairs)

    # Forming C-hetero bonds and breaking C-H bonds are oxidation-like.
    score += float(formed_het)
    score += float(broken_ch)
    # Breaking C-hetero bonds and forming C-H bonds are reduction-like.
    score -= float(broken_het)
    score -= float(formed_ch)

    if formed_het:
        signals.append(f"formed_c_hetero={formed_het}")
    if broken_het:
        signals.append(f"broken_c_hetero={broken_het}")
    if formed_ch:
        signals.append("formed_c_h")
    if broken_ch:
        signals.append("broken_c_h")

    has_evidence = bool(formed_pairs or broken_pairs)
    if not has_evidence:
        return {
            "classification": "uncertain",
            "confidence": 0.2,
            "score": 0.0,
            "reasons": ["missing_bond_change_evidence"],
            "signals": [],
        }

    if abs(score) < 0.5:
        classification = "redox_neutral"
    elif score > 0:
        classification = "net_oxidation"
    else:
        classification = "net_reduction"

    # Confidence is higher when there are explicit C-H changes or multiple
    # oxidation-state-relevant bond changes.
    informative = formed_het + broken_het + formed_ch + broken_ch
    confidence = 0.55 + min(0.35, 0.1 * max(0, informative - 1))
    confidence = round(max(0.2, min(0.95, confidence)), 2)

    reasons: List[str] = []
    if classification == "redox_neutral":
        reasons.append("oxidation_and_reduction_signals_balanced_or_absent")
    elif classification == "net_oxidation":
        reasons.append("oxidation_signals_exceed_reduction_signals")
    else:
        reasons.append("reduction_signals_exceed_oxidation_signals")

    return {
        "classification": classification,
        "confidence": confidence,
        "score": round(score, 2),
        "reasons": reasons,
        "signals": signals,
    }


def _count_molecules_matching_smarts(
    smiles_list: Iterable[str],
    smarts_list: Iterable[str],
) -> int:
    if not rdkit_helpers.rdkit_available():
        return 0
    try:
        from chemtools.util.smarts_cache import compile_smarts
    except Exception:
        return 0
    patterns = []
    for smarts in smarts_list:
        patt = compile_smarts(smarts, validate=False)
        if patt is not None:
            patterns.append(patt)
    if not patterns:
        return 0

    hits = 0
    for smiles in smiles_list or []:
        mol = rdkit_helpers.parse_smiles(smiles)
        if mol is None:
            continue
        try:
            if any(mol.HasSubstructMatch(p) for p in patterns):
                hits += 1
        except Exception:
            continue
    return hits


def _infer_event_families(events: Iterable[Dict[str, Any]]) -> Set[str]:
    kinds = {
        str(ev.get("kind")).strip()
        for ev in (events or [])
        if isinstance(ev, dict) and str(ev.get("kind") or "").strip()
    }
    families: Set[str] = set()
    has_benzyl_o_alkylation = "benzyl_o_alkylation_like" in kinds
    if has_benzyl_o_alkylation:
        families.add("o_alkylation")

    if "ester_hydrolysis_like" in kinds:
        families.add("hydrolysis")
    if "amidation_like" in kinds:
        families.add("acyl_transfer")
    if "ring_closure_or_annulation" in kinds:
        families.add("annulation")

    bond_formation_kinds = {
        "c_n_bond_formation",
        "c_o_bond_formation",
        "c_s_bond_formation",
        "c_c_bond_formation",
    }
    has_bond_formation = bool(kinds & bond_formation_kinds)
    has_displacement = "leaving_group_displacement" in kinds
    if not has_benzyl_o_alkylation:
        if has_displacement and has_bond_formation:
            families.add("substitution")
        elif has_displacement:
            families.add("displacement")
        elif has_bond_formation:
            families.add("bond_formation")
    return families


def format_multi_event_signature(events: Iterable[Dict[str, Any]]) -> str:
    kinds = {
        str(ev.get("kind")).strip()
        for ev in (events or [])
        if isinstance(ev, dict) and str(ev.get("kind") or "").strip()
    }
    ordered_codes: List[str] = []
    for kind in _EVENT_SIGNATURE_PRIORITY:
        if kind not in kinds:
            continue
        code = _EVENT_SIGNATURE_CODE.get(kind)
        if code and code not in ordered_codes:
            ordered_codes.append(code)
    if len(ordered_codes) < 2:
        return ""
    return "+".join(ordered_codes[:4])


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

    reactant_ester_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["ester"]],
    )
    product_ester_count = _count_molecules_matching_smarts(
        product_smiles,
        [_PATTERN_SMARTS["ester"]],
    )
    reactant_acid_like_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["carboxy_acid"], _PATTERN_SMARTS["carboxylate"]],
    )
    product_acid_like_count = _count_molecules_matching_smarts(
        product_smiles,
        [_PATTERN_SMARTS["carboxy_acid"], _PATTERN_SMARTS["carboxylate"]],
    )
    if (
        reactant_ester_count > product_ester_count
        and product_acid_like_count > reactant_acid_like_count
    ):
        events.append(
            _event(
                "ester_hydrolysis_like",
                0.82,
                {
                    "reactant_ester_count": reactant_ester_count,
                    "product_ester_count": product_ester_count,
                    "reactant_acid_like_count": reactant_acid_like_count,
                    "product_acid_like_count": product_acid_like_count,
                },
            )
        )

    reactant_phenol_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["phenol"]],
    )
    reactant_benzyl_halide_count = _count_molecules_matching_smarts(
        reactant_smiles,
        [_PATTERN_SMARTS["benzyl_halide"]],
    )
    product_benzyl_aryl_ether_count = _count_molecules_matching_smarts(
        product_smiles,
        [_PATTERN_SMARTS["benzyl_aryl_ether"]],
    )
    if (
        reactant_phenol_count > 0
        and reactant_benzyl_halide_count > 0
        and product_benzyl_aryl_ether_count > 0
    ):
        events.append(
            _event(
                "benzyl_o_alkylation_like",
                0.84,
                {
                    "reactant_phenol_count": reactant_phenol_count,
                    "reactant_benzyl_halide_count": reactant_benzyl_halide_count,
                    "product_benzyl_aryl_ether_count": product_benzyl_aryl_ether_count,
                },
            )
        )

    event_families = _infer_event_families(events)
    has_complex_family = bool(event_families & {"hydrolysis", "acyl_transfer", "annulation"})
    if len(event_families) >= 3 or (len(event_families) >= 2 and has_complex_family):
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
    redox_assessment = _compute_redox_assessment(formed_pairs, broken_pairs)

    return {
        "events": events,
        "anomalies": anomalies,
        "event_families": sorted(event_families),
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
        "redox_assessment": redox_assessment,
    }


__all__ = [
    "format_multi_event_signature",
    "summarize_reaction_events",
]
