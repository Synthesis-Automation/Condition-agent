"""
Reaction-level feature formatting and bundling.

Handles reaction type detection, reactant/product processing, and reaction key generation.
"""

from __future__ import annotations

from functools import lru_cache
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.util import rdkit_helpers
from chemtools.smiles import normalize_reaction

from ..analysis.reaction_context import classify_reactants_with_context, get_reactant_summary
from ..analysis.feasibility import analyze_snar_feasibility
from ..analysis.reaction_record import ReactionRecord

from .molecule import build_molecule_bundle, to_bool
from .aggregation import aggregate_reaction_features, infer_intramolecular
from .utils import extract_motif_ids, normalize_motif_id
from .simplified import build_core_reaction, build_extended_reaction
from .reaction_events import summarize_reaction_events
from ..spectator_rank import rank_spectator_groups


_UNCLASSIFIED_REACTANT_MOTIF = "Unclassified-Reactant"
_HYBRID_MAPPING_MIN_CONFIDENCE = 0.65
_LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD = 0.60


def _to_float_or_default(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def _get_reaction_type_id_and_confidence(
    reaction_type: Any,
) -> Tuple[Optional[str], float]:
    if isinstance(reaction_type, dict):
        rid = reaction_type.get("reaction_type")
        conf = _to_float_or_default(reaction_type.get("confidence"), 0.0)
    else:
        rid = reaction_type
        conf = 0.0
    rid_text = str(rid or "").strip()
    if not rid_text or rid_text.lower() == "unknown":
        return None, conf
    return rid_text, conf


def _set_reaction_type_payload(
    reaction_type: Any,
    reaction_type_id: str,
    confidence: float,
) -> Dict[str, Any]:
    payload: Dict[str, Any] = (
        dict(reaction_type) if isinstance(reaction_type, dict) else {}
    )
    payload["reaction_type"] = reaction_type_id
    payload["name"] = reaction_type_id
    payload["confidence"] = max(0.0, min(1.0, _to_float_or_default(confidence, 0.0)))
    return payload


def _normalize_llm_assist_options(options: Dict[str, Any]) -> Dict[str, Any]:
    cfg = options.get("llm_assist")
    defaults = {
        "enabled": False,
        "provider": None,
        "model": None,
        "temperature": 0.0,
        "max_tokens": 700,
        "timeout": 60,
        "only_on_uncertain": True,
        "confidence_threshold": _LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD,
        "require_crk_validation": True,
    }
    if isinstance(cfg, bool):
        defaults["enabled"] = cfg
        return defaults
    if not isinstance(cfg, dict):
        return defaults
    defaults["enabled"] = to_bool(cfg.get("enabled"), default=True)
    defaults["provider"] = (
        str(cfg.get("provider")).strip() if cfg.get("provider") is not None else None
    )
    defaults["model"] = (
        str(cfg.get("model")).strip() if cfg.get("model") is not None else None
    )
    defaults["temperature"] = _to_float_or_default(cfg.get("temperature"), 0.0)
    defaults["max_tokens"] = int(cfg.get("max_tokens") or 700)
    defaults["timeout"] = int(cfg.get("timeout") or 60)
    defaults["only_on_uncertain"] = to_bool(cfg.get("only_on_uncertain"), default=True)
    defaults["confidence_threshold"] = _to_float_or_default(
        cfg.get("confidence_threshold"),
        _LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD,
    )
    defaults["require_crk_validation"] = to_bool(
        cfg.get("require_crk_validation"),
        default=True,
    )
    return defaults


def _is_reaction_uncertain_for_llm_assist(
    reaction_type: Any,
    detection_payload: Dict[str, Any],
    reaction_key: Optional[str],
    confidence_threshold: float,
) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    rt_id, confidence = _get_reaction_type_id_and_confidence(reaction_type)
    if not rt_id:
        reasons.append("unknown_reaction_type")
    if confidence < max(0.0, confidence_threshold):
        reasons.append("low_confidence")
    if isinstance(detection_payload.get("mapping_warning"), dict):
        reasons.append("mapping_warning")
    if not str(reaction_key or "").strip():
        reasons.append("missing_reaction_key")
    key_quality = detection_payload.get("reaction_key_quality")
    if isinstance(key_quality, dict):
        level = str(key_quality.get("level") or "").strip().lower()
        score = _to_float_or_default(key_quality.get("score_0_1"), 1.0)
        if level == "low" or score < 0.45:
            reasons.append("low_reaction_key_quality")
        elif level == "medium":
            reasons.append("medium_reaction_key_quality")
    validation_payload = detection_payload.get("validation")
    if isinstance(validation_payload, dict):
        validated = str(validation_payload.get("validated_detection") or "").strip()
        if not validated or validated.lower() == "unknown":
            reasons.append("unknown_validated_detection")
    return bool(reasons), reasons


def _run_primary_reaction_type_detection(
    reaction_smiles: str,
) -> Dict[str, Any]:
    """
    Run the canonical detection API and return normalized payload.

    This is optional because it can be more expensive than local validation.
    """
    try:
        from chemtools.detection import detect_reaction_type
    except Exception as exc:
        return {
            "status": "error",
            "error": f"detection_import_error: {exc}",
            "summary": None,
            "matches": [],
        }
    try:
        result = detect_reaction_type(reaction_smiles)
    except Exception as exc:
        return {
            "status": "error",
            "error": f"detection_runtime_error: {exc}",
            "summary": None,
            "matches": [],
        }

    summary = format_reaction_type_summary(result)
    matches_payload: List[Dict[str, Any]] = []
    for match in getattr(result, "matches", []) or []:
        slot_evidence: Dict[str, Any] = {}
        electrophile = list(getattr(match, "electrophile", []) or [])
        nucleophile = list(getattr(match, "nucleophile", []) or [])
        product = list(getattr(match, "product", []) or [])
        if electrophile:
            slot_evidence["electrophile"] = electrophile
        if nucleophile:
            slot_evidence["nucleophile"] = nucleophile
        if product:
            slot_evidence["product"] = product
        matches_payload.append(
            {
                "reaction_type": str(getattr(match, "reaction_type", "") or ""),
                "name": str(getattr(match, "reaction_type", "") or ""),
                "confidence": _to_float_or_default(getattr(match, "confidence", 0.0), 0.0),
                "slot_evidence": slot_evidence,
            }
        )
    return {
        "status": "ok",
        "error": getattr(result, "error", None),
        "summary": summary,
        "matches": matches_payload,
    }


def _run_llm_reaction_assist(
    context: Dict[str, Any],
    llm_assist: Dict[str, Any],
) -> Dict[str, Any]:
    try:
        from llmtools.reaction_featurization_review import (
            LLMReactionFeaturizationOptions,
            review_reaction_featurization,
        )
    except Exception as exc:
        return {
            "enabled": True,
            "status": "error",
            "error": f"llmtools unavailable: {exc}",
        }

    opts = LLMReactionFeaturizationOptions(
        enabled=True,
        provider=llm_assist.get("provider"),
        model=llm_assist.get("model"),
        temperature=_to_float_or_default(llm_assist.get("temperature"), 0.0),
        max_tokens=int(llm_assist.get("max_tokens") or 700),
        timeout=int(llm_assist.get("timeout") or 60),
    )
    return review_reaction_featurization(context, opts)


def _validate_llm_suggested_type_with_crk(
    suggested_reaction_type: str,
    reaction_key_raw: Optional[str],
) -> Tuple[bool, str]:
    if not str(reaction_key_raw or "").strip():
        return True, "no_reaction_key_for_validation"
    try:
        from .detection_validation import validate_detection_with_crk_key
    except Exception as exc:
        return False, f"validation_import_error: {exc}"
    try:
        validated = validate_detection_with_crk_key(
            initial_detection=suggested_reaction_type,
            initial_confidence=0.9,
            reaction_key=str(reaction_key_raw),
            include_evidence=False,
        )
    except Exception as exc:
        return False, f"validation_error: {exc}"
    validated_type = _resolve_reaction_type_id(validated.get("reaction_type"))
    if validated_type == suggested_reaction_type:
        return True, "validated_with_crk"
    return False, f"validation_mismatch:{validated_type or 'Unknown'}"


def _hybrid_mapping_is_unreliable(analysis: Dict[str, Any]) -> bool:
    """
    Flag low-confidence hybrid mapping when RXNMapper and MCS disagree.

    In this case, mapped bond deltas are often noisy and can remove true
    reaction-center signals downstream after bond sanitization.
    """
    if not isinstance(analysis, dict):
        return False
    if str(analysis.get("method") or "").strip() != "hybrid":
        return False
    agreement = (analysis.get("agreement") or {}).get("rxnmapper_vs_mcs")
    if agreement is not False:
        return False
    combined_conf = _to_float_or_default(analysis.get("combined_confidence"), 0.0)
    return combined_conf < _HYBRID_MAPPING_MIN_CONFIDENCE


def _preferred_bond_change_result(analysis: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    if not isinstance(analysis, dict):
        return None
    for key in ("recommended_result", "rxnmapper_result", "manual_result"):
        candidate = analysis.get(key)
        if isinstance(candidate, dict) and candidate.get("success"):
            return candidate
    return None


def _get_bond_change_analysis_with_quality(
    reaction_smiles: str,
) -> Tuple[Optional[Dict[str, Any]], bool, Optional[Dict[str, Any]]]:
    """
    Return preferred bond-change analysis with reliability metadata.

    Returns:
        (preferred_result, unreliable, raw_hybrid_analysis)
    """
    if not reaction_smiles or not rdkit_helpers.rdkit_available():
        return None, False, None
    try:
        from chemtools._atom_mapping import analyze_bond_changes_hybrid
    except Exception:
        return None, False, None
    try:
        analysis = analyze_bond_changes_hybrid(reaction_smiles, use_mcs=True)
    except Exception:
        return None, False, None
    if not analysis or not analysis.get("success"):
        return None, False, analysis if isinstance(analysis, dict) else None
    preferred = _preferred_bond_change_result(analysis)
    if not preferred:
        return None, False, analysis
    return preferred, _hybrid_mapping_is_unreliable(analysis), analysis


def get_crk_options() -> Dict[str, Any]:
    """Return standardized options for CRK-v1 reaction key generation."""
    return {
        "include_roles": False,
        "include_agent_roles": False,
        "include_reaction_type": False,
        "motif_site_filter": "substituent",
        "confirm_coupling_products": True,
        "discovery_mode": False,
        "reactant_coverage_guard": True,
    }


def _bundle_has_taxonomy_motif(bundle: Dict[str, Any]) -> bool:
    for key in ("motifs", "context_motifs"):
        for motif in bundle.get(key) or []:
            cid = ""
            if isinstance(motif, dict):
                cid = str(motif.get("compound_id") or motif.get("id") or "").strip()
            else:
                cid = str(motif).strip()
            if not cid or cid == "unknown":
                continue
            if cid.startswith("Unclassified-"):
                continue
            return True
    return False


def _ensure_reactant_coverage(
    reactant_bundles: List[Dict[str, Any]],
    *,
    enabled: bool,
) -> None:
    if not enabled:
        return
    for bundle in reactant_bundles:
        if _bundle_has_taxonomy_motif(bundle):
            continue
        motifs = bundle.get("motifs")
        if not isinstance(motifs, list):
            motifs = []
            bundle["motifs"] = motifs
        motifs.append(
            {
                "compound_id": _UNCLASSIFIED_REACTANT_MOTIF,
                "undocumented": True,
                "coverage_guard": True,
                "reason": "no_taxonomy_motif_detected",
            }
        )
        ranked = bundle.get("ranked_motifs")
        if isinstance(ranked, list):
            if _UNCLASSIFIED_REACTANT_MOTIF not in ranked:
                ranked.append(_UNCLASSIFIED_REACTANT_MOTIF)


def format_reaction_type_summary(detection: Any) -> Dict[str, Any]:
    """Extract reaction type information with alternatives."""
    matches = detection.matches if detection else []
    if not matches:
        return {"reaction_type": "Unknown", "confidence": 0.0, "slot_evidence": {}}
    
    best = matches[0]
    # Construct slot_evidence from new API structure
    slot_evidence = {}
    if best.electrophile:
        slot_evidence["electrophile"] = best.electrophile
    if best.nucleophile:
        slot_evidence["nucleophile"] = best.nucleophile
    if best.product:
        slot_evidence["product"] = best.product
    
    result = {
        "reaction_type": best.reaction_type,
        "name": best.reaction_type,  # Use reaction_type as name for compatibility
        "category": "coupling",  # Default category - could be improved
        "confidence": best.confidence,
        "slot_evidence": slot_evidence,
    }
    
    # Add alternatives if available (top 3 total)
    if len(matches) > 1:
        alts = []
        for alt in matches[1:3]:
            alts.append({
                "reaction_type": alt.reaction_type,
                "name": alt.reaction_type,
                "confidence": alt.confidence
            })
        result["alternatives"] = alts
        
    return result


def _map_atom_labels(mapped_smiles: str) -> Dict[int, str]:
    """Build a map_num -> label mapping from a mapped reaction SMILES."""
    if not mapped_smiles or ">>" not in mapped_smiles:
        return {}
    labels: Dict[int, str] = {}
    reactants_smiles, products_smiles = mapped_smiles.split(">>", 1)
    for side in (reactants_smiles, products_smiles):
        for part in side.split("."):
            mol = rdkit_helpers.parse_smiles(part)
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if not map_num or map_num in labels:
                    continue
                symbol = atom.GetSymbol()
                if atom.GetIsAromatic():
                    labels[map_num] = f"{symbol}(ar)"
                else:
                    labels[map_num] = symbol
    return labels


def _format_bond_tokens(
    bonds: Iterable[Any],
    labels: Dict[int, str],
) -> List[str]:
    tokens: List[str] = []
    seen: Set[Tuple[Any, ...]] = set()
    for bond in bonds or []:
        if not isinstance(bond, (tuple, list)) or len(bond) < 2:
            continue
        a, b = bond[0], bond[1]
        if isinstance(a, int):
            a_label = labels.get(a, f"#{a}")
        else:
            a_label = str(a).split()[0]
        if isinstance(b, int):
            b_label = labels.get(b, f"#{b}")
        else:
            b_label = str(b).split()[0]
        pair = tuple(sorted((a_label, b_label)))
        if isinstance(a, int) and isinstance(b, int):
            bond_id = ("map", min(a, b), max(a, b))
        else:
            bond_id = ("label", pair[0], pair[1])
        if bond_id in seen:
            continue
        seen.add(bond_id)
        tokens.append(f"{pair[0]}-{pair[1]}")
    return tokens


def _get_bond_change_analysis(
    reaction_smiles: str,
    *,
    allow_unreliable: bool = False,
) -> Optional[Dict[str, Any]]:
    result, unreliable, _raw = _get_bond_change_analysis_with_quality(reaction_smiles)
    if unreliable and not allow_unreliable:
        # Prefer no bond key over low-confidence mapped bond deltas.
        # This preserves motif-based CRK reactants when mapping is noisy.
        return None
    return result


@lru_cache(maxsize=1)
def _load_compound_logic_sets() -> Dict[str, Set[str]]:
    path = Path(__file__).resolve().parent.parent.parent / "taxonomy" / "data" / "compound_logic.json"
    if not path.exists():
        return {}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            payload = json.load(handle)
    except Exception:
        return {}
    motif_sets: Dict[str, Set[str]] = {}
    for set_name, set_data in (payload.get("motif_sets", {}) or {}).items():
        members = set_data.get("members", []) or []
        motif_sets[set_name] = set(str(m) for m in members if m)
    return motif_sets


@lru_cache(maxsize=1)
def _load_group_element_map() -> Dict[str, Set[str]]:
    try:
        from chemtools.taxonomy import loader as taxonomy_loader
    except Exception:
        return {}
    payload = taxonomy_loader.load_group_logic()
    if not payload:
        return {}
    mapping = payload.get("group_elements", {}) or {}
    group_element_map: Dict[str, Set[str]] = {}
    for group_id, elements in mapping.items():
        if not group_id or not elements:
            continue
        if isinstance(elements, (list, tuple, set)):
            values = {str(el) for el in elements if el}
        else:
            values = {str(elements)}
        if values:
            group_element_map[str(group_id)] = values
    return group_element_map


@lru_cache(maxsize=1)
def _load_group_inference_rules() -> Dict[str, Any]:
    try:
        from chemtools.taxonomy import loader as taxonomy_loader
    except Exception:
        return {}
    payload = taxonomy_loader.load_group_logic()
    if not payload:
        return {}
    rules = payload.get("product_inference", {}) or {}
    if not isinstance(rules, dict):
        return {}
    return rules


@lru_cache(maxsize=1)
def _load_reaction_catalog_data() -> Tuple[Dict[str, Any], Dict[str, str]]:
    try:
        from chemtools.taxonomy.reaction_catalog import load_reaction_catalog
    except Exception:
        return {}, {}
    try:
        definitions, alias_map = load_reaction_catalog()
    except Exception:
        return {}, {}
    return definitions, alias_map


def _resolve_reaction_type_id(reaction_type: Optional[str]) -> Optional[str]:
    if not reaction_type:
        return None
    definitions, alias_map = _load_reaction_catalog_data()
    if not definitions:
        return None
    label = str(reaction_type).strip()
    if not label:
        return None
    if label in definitions:
        return label
    resolved = alias_map.get(label.lower())
    if resolved and resolved in definitions:
        return resolved
    return None


def _taxonomy_allowed_product_motifs(reaction_type: Optional[str]) -> Set[str]:
    resolved = _resolve_reaction_type_id(reaction_type)
    if not resolved:
        return set()
    definitions, _ = _load_reaction_catalog_data()
    definition = definitions.get(resolved)
    if not definition:
        return set()
    allowed: Set[str] = set()
    for slot_req in definition.products.values():
        for motif in slot_req.allowed:
            motif_id = str(motif).strip()
            if motif_id:
                allowed.add(motif_id)
    return allowed


def format_bond_change_key(
    reaction_smiles: str,
    *,
    analysis: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Generate a bond-change key from atom mapping (POC)."""
    result = analysis or _get_bond_change_analysis(reaction_smiles)
    if not result:
        return None
    mapped_smiles = result.get("mapped_smiles") or ""
    labels = _map_atom_labels(mapped_smiles)
    broken = _format_bond_tokens(result.get("broken_bonds"), labels)
    formed = _format_bond_tokens(result.get("formed_bonds"), labels)
    if broken and formed:
        # Remove any bonds that appear as both broken and formed (mapping artifacts).
        overlap = set(broken) & set(formed)
        if overlap:
            from collections import Counter
            broken_counts = Counter(broken)
            formed_counts = Counter(formed)
            for token in overlap:
                drop = min(broken_counts.get(token, 0), formed_counts.get(token, 0))
                if drop <= 0:
                    continue
                broken_counts[token] -= drop
                formed_counts[token] -= drop
            broken = [t for t, count in broken_counts.items() for _ in range(count) if count > 0]
            formed = [t for t, count in formed_counts.items() for _ in range(count) if count > 0]
    if not broken and not formed:
        return None
    parts: List[str] = []
    if broken:
        parts.append("break: " + "; ".join(broken))
    if formed:
        parts.append("form: " + "; ".join(formed))
    return " | ".join(parts)


def _extract_bond_section(bond_key: Optional[str], *, section: str) -> List[str]:
    if not bond_key:
        return []
    marker = f"{section}: "
    parts = [p.strip() for p in bond_key.split(" | ") if p.strip()]
    for part in parts:
        if part.startswith(marker):
            payload = part[len(marker):].strip()
            if not payload:
                return []
            return [tok.strip() for tok in payload.split(";") if tok.strip()]
    return []


def _is_carbon_label(label: str) -> bool:
    return label == "C" or label.startswith("C(")


def _is_sulfur_label(label: str) -> bool:
    return label == "S" or label.startswith("S(")


def _parse_element(label: str) -> str:
    label = label.strip()
    if label.startswith("#"):
        return ""
    if "(" in label:
        label = label.split("(", 1)[0]
    return label


def _is_aromatic_label(label: str) -> bool:
    return label.endswith("(ar)")


def _bond_elements_from_key(bond_key: Optional[str]) -> Set[str]:
    elements: Set[str] = set()
    for section in ("break", "form"):
        for token in _extract_bond_section(bond_key, section=section):
            if "-" not in token:
                continue
            left, right = [t.strip() for t in token.split("-", 1)]
            left_el = _parse_element(left)
            right_el = _parse_element(right)
            if left_el:
                elements.add(left_el)
            if right_el:
                elements.add(right_el)
    return elements


def _label_formed_bonds(
    formed: Iterable[str],
    product_motifs_reactive: Iterable[str],
) -> List[str]:
    """Attach product motif labels to formed bond tokens when possible."""
    tokens = [str(t) for t in formed if t]
    motifs = [str(m) for m in product_motifs_reactive if m]
    if not tokens or not motifs:
        return []

    cc_labels = [m for m in motifs if m in {"Ar-Ar", "Ar-Alkyl", "Ar-Alkenyl", "Ar-Alkynyl"}]
    cn_labels = [m for m in motifs if m.startswith("Ar-N") or m == "Ar-AromN"]
    co_labels = [m for m in motifs if m == "Ar-OR"]
    cs_labels = [m for m in motifs if m == "Ar-SR"]

    cc_idx = 0
    cn_idx = 0
    co_idx = 0
    cs_idx = 0

    labeled: List[str] = []
    for token in tokens:
        if "-" not in token:
            labeled.append(token)
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        left_ar = _is_aromatic_label(left)
        right_ar = _is_aromatic_label(right)
        elements = {left_el, right_el}
        is_aryl_carbon = (left_el == "C" and left_ar) or (right_el == "C" and right_ar)

        label = None
        if elements == {"C"} and is_aryl_carbon:
            if cc_idx < len(cc_labels):
                label = cc_labels[cc_idx]
                cc_idx += 1
        elif "N" in elements and "C" in elements and is_aryl_carbon:
            if cn_idx < len(cn_labels):
                label = cn_labels[cn_idx]
                cn_idx += 1
        elif "O" in elements and "C" in elements and is_aryl_carbon:
            if co_idx < len(co_labels):
                label = co_labels[co_idx]
                co_idx += 1
        elif "S" in elements and "C" in elements and is_aryl_carbon:
            if cs_idx < len(cs_labels):
                label = cs_labels[cs_idx]
                cs_idx += 1

        if label:
            labeled.append(f"{token}[{label}]")
        else:
            labeled.append(token)
    return labeled


@lru_cache(maxsize=1)
def _load_group_sets() -> Dict[str, Any]:
    try:
        from chemtools.taxonomy import loader as taxonomy_loader
    except Exception:
        return {}
    payload = taxonomy_loader.load_group_logic()
    if not payload:
        return {}
    group_sets = payload.get("group_sets", {}) or {}
    return group_sets if isinstance(group_sets, dict) else {}


def _expand_group_set(
    set_id: str,
    group_sets: Dict[str, Any],
    *,
    _seen: Optional[Set[str]] = None,
) -> Set[str]:
    if not group_sets or not set_id:
        return set()
    if _seen is None:
        _seen = set()
    if set_id in _seen:
        return set()
    _seen.add(set_id)
    set_data = group_sets.get(set_id, {}) or {}
    members = set_data.get("members", []) or []
    expanded: Set[str] = set()
    for member in members:
        if not member:
            continue
        member_str = str(member)
        if member_str in group_sets:
            expanded.update(_expand_group_set(member_str, group_sets, _seen=_seen))
        else:
            expanded.add(member_str)
    return expanded


@lru_cache(maxsize=1)
def _load_leaving_groups() -> Set[str]:
    group_sets = _load_group_sets()
    if not group_sets:
        return set()
    return _expand_group_set("LeavingGroup", group_sets)


def _reacted_has_leaving_groups(reacted: Iterable[str]) -> bool:
    leaving_groups = _load_leaving_groups()
    if not leaving_groups:
        return False
    for motif in reacted:
        motif_str = str(motif)
        for group_id in leaving_groups:
            if motif_str.endswith(group_id):
                return True
    return False


def _nucleophile_elements_from_motifs(
    reacted: Iterable[str],
    group_element_map: Dict[str, Set[str]],
) -> Set[str]:
    elements: Set[str] = set()
    for motif in reacted:
        group_id = _match_group_id(str(motif), group_element_map)
        if not group_id:
            continue
        elements.update(group_element_map.get(group_id, set()))
    elements.discard("C")
    elements.discard("H")
    return elements


def _sanitize_bond_key(
    bond_key: Optional[str],
    reacted: Iterable[str],
    *,
    group_element_map: Dict[str, Set[str]],
) -> Optional[str]:
    if not bond_key:
        return None
    broken = _extract_bond_section(bond_key, section="break")
    formed = _extract_bond_section(bond_key, section="form")
    if not broken and not formed:
        return None
    if not broken and formed and _reacted_has_leaving_groups(reacted):
        nucleophile_elements = _nucleophile_elements_from_motifs(reacted, group_element_map)
        if nucleophile_elements:
            filtered = []
            for token in formed:
                if "-" not in token:
                    continue
                left, right = [t.strip() for t in token.split("-", 1)]
                left_el = _parse_element(left)
                right_el = _parse_element(right)
                if left_el in nucleophile_elements or right_el in nucleophile_elements:
                    filtered.append(token)
            formed = filtered
        else:
            formed = []
        if not formed:
            return None
    parts: List[str] = []
    if broken:
        parts.append("break: " + "; ".join(broken))
    if formed:
        parts.append("form: " + "; ".join(formed))
    return " | ".join(parts) if parts else None


def _match_group_id(motif_id: str, group_element_map: Dict[str, Set[str]]) -> Optional[str]:
    motif_id = str(motif_id)
    for group_id in sorted(group_element_map.keys(), key=len, reverse=True):
        if motif_id.endswith(group_id):
            return group_id
    return None


def _filter_reactants_for_crk(
    reacted: Iterable[str],
    bond_key: Optional[str],
    spectators: Optional[Iterable[str]] = None,
) -> List[str]:
    scaffold_ids: Set[str] = set()
    try:
        from .aggregation import load_scaffold_motif_ids
        scaffold_ids = load_scaffold_motif_ids()
    except Exception:
        scaffold_ids = set()

    if not bond_key:
        return sorted(str(r) for r in reacted if r and str(r) not in scaffold_ids)
    elements = _bond_elements_from_key(bond_key)
    if not elements:
        return sorted(str(r) for r in reacted if r and str(r) not in scaffold_ids)
    elements_non_c = {el for el in elements if el != "C"}

    # If N/O/S bonds are formed, make sure N/O/S nucleophiles are not dropped
    # just because they have -H suffixes (e.g., AromN-H).
    formed_tokens = _extract_bond_section(bond_key, section="form")
    formed_elements: Set[str] = set()
    for token in formed_tokens:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        if left_el:
            formed_elements.add(left_el)
        if right_el:
            formed_elements.add(right_el)

    bond_tokens = _extract_bond_section(bond_key, section="form") + _extract_bond_section(
        bond_key, section="break"
    )

    def _parse_bond_token(token: str) -> Optional[Tuple[str, bool, str, bool]]:
        if "-" not in token:
            return None
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        if not left_el or not right_el:
            return None
        return (left_el, _is_aromatic_label(left), right_el, _is_aromatic_label(right))

    parsed_bonds = [entry for entry in (_parse_bond_token(t) for t in bond_tokens) if entry]

    def _has_c_c_bond() -> bool:
        for left_el, _left_ar, right_el, _right_ar in parsed_bonds:
            if left_el == "C" and right_el == "C":
                return True
        return False

    def _has_c_c_cleavage() -> bool:
        for token in _extract_bond_section(bond_key, section="break"):
            if "-" not in token:
                continue
            left, right = [t.strip() for t in token.split("-", 1)]
            if _parse_element(left) == "C" and _parse_element(right) == "C":
                return True
        return False

    def _has_aryl_c_c_formation() -> bool:
        for token in _extract_bond_section(bond_key, section="form"):
            if "-" not in token:
                continue
            left, right = [t.strip() for t in token.split("-", 1)]
            left_el = _parse_element(left)
            right_el = _parse_element(right)
            left_ar = _is_aromatic_label(left)
            right_ar = _is_aromatic_label(right)
            if left_el == "C" and right_el == "C" and (left_ar or right_ar):
                return True
        return False

    def _has_non_aromatic_c_bond(other_elements: Set[str]) -> bool:
        for left_el, left_ar, right_el, right_ar in parsed_bonds:
            if left_el == "C" and not left_ar and right_el in other_elements:
                return True
            if right_el == "C" and not right_ar and left_el in other_elements:
                return True
        return False

    group_element_map = _load_group_element_map()

    reacted_ids = {str(r) for r in reacted if r}
    spectator_ids = {str(s) for s in (spectators or []) if s}

    # Ensure aromatic C-H partner for Minisci-like decarboxylative C-C(ar)
    # formation so CRK reactants stay mechanistic even if C-H motifs are
    # simultaneously marked as spectators by motif-delta logic.
    promoted_from_spectators: List[str] = []
    must_keep_aromatic_ch: Set[str] = set()
    logic_sets = _load_compound_logic_sets()
    carboxylic_acids = logic_sets.get("carboxylic_acids", set())
    has_decarboxylative_partner = any(
        (m in carboxylic_acids) or m.endswith("-CO2H") or m.endswith("-COOH")
        for m in reacted_ids
    )
    has_aryl_c_c_formation = False
    for token in formed_tokens:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        left_ar = _is_aromatic_label(left)
        right_ar = _is_aromatic_label(right)
        if left_el == "C" and right_el == "C" and (left_ar or right_ar):
            has_aryl_c_c_formation = True
            break
    if has_decarboxylative_partner and has_aryl_c_c_formation:
        for candidate in ("HeteroAr-H", "AromN-H", "Ar-H"):
            if candidate in reacted_ids or candidate in spectator_ids:
                must_keep_aromatic_ch.add(candidate)
                if candidate in spectator_ids and candidate not in reacted_ids:
                    promoted_from_spectators.append(candidate)
                break

    filtered: List[str] = []
    for motif in reacted:
        if not motif:
            continue
        if str(motif) in scaffold_ids:
            continue
        motif_str = str(motif)
        scaffold = motif_str.split("-", 1)[0] if "-" in motif_str else motif_str
        if motif_str in must_keep_aromatic_ch:
            filtered.append(motif_str)
            continue
        # Keep N/O/S nucleophiles when corresponding bonds form.
        keep_by_formed = False
        if "N" in formed_elements and (
            motif_str.endswith("-NH2")
            or motif_str.endswith("-NHR")
            or motif_str.endswith("-NR2")
            or motif_str.endswith("-NHCOR")
            or motif_str == "AromN-H"
        ):
            keep_by_formed = True
        if "O" in formed_elements and (motif_str.endswith("-OH") or motif_str.endswith("-OR")):
            keep_by_formed = True
        if "S" in formed_elements and (motif_str.endswith("-SH") or motif_str.endswith("-SR")):
            keep_by_formed = True
        if scaffold.startswith("Alkenyl") or scaffold.startswith("Alkynyl"):
            if _has_c_c_bond():
                keep_by_formed = True

        group_id = _match_group_id(motif_str, group_element_map)
        keep_by_elements = False
        if group_id:
            group_elements = group_element_map.get(group_id)
            if group_elements:
                if elements_non_c:
                    if elements_non_c.intersection(group_elements):
                        keep_by_elements = True
                elif elements.intersection(group_elements):
                    keep_by_elements = True

        # For carbonyl-derived groups, element overlap alone is too permissive.
        # Require explicit nucleophile logic to keep these when they are spectators.
        if group_id in {
            "-CO2R",
            "-CO2H",
            "-COR",
            "-CONH2",
            "-CONHR",
            "-CONR2",
        }:
            # Carbonyl-derived motifs should only be kept if a non-aromatic carbonyl bond changes.
            if (
                group_id == "-CO2H"
                and _has_c_c_cleavage()
                and _has_aryl_c_c_formation()
            ):
                keep_by_elements = True
            elif not _has_non_aromatic_c_bond({"O", "N", "Cl", "Br", "I", "F", "S"}):
                keep_by_elements = False

        # Drop pure spectators unless they align with bond elements or formed nucleophiles.
        if motif_str in spectator_ids and not (keep_by_formed or keep_by_elements):
            continue

        if keep_by_formed:
            filtered.append(motif_str)
            continue
        if not group_id:
            filtered.append(motif_str)
            continue
        if keep_by_elements:
            filtered.append(motif_str)
    filtered_sorted = sorted(filtered)
    if filtered_sorted:
        if promoted_from_spectators:
            filtered_sorted = sorted(set(filtered_sorted) | set(promoted_from_spectators))
        return filtered_sorted
    fallback = [str(r) for r in reacted if r and str(r) not in scaffold_ids]
    if fallback:
        fallback_sorted = sorted(fallback)
        if promoted_from_spectators:
            fallback_sorted = sorted(set(fallback_sorted) | set(promoted_from_spectators))
        return fallback_sorted
    if promoted_from_spectators:
        return sorted(set(promoted_from_spectators))
    return sorted(str(r) for r in reacted if r)


def _strip_atom_mapping(smiles: str) -> str:
    return re.sub(r":\d+", "", smiles)


def _build_map_info(
    mapped_smiles: str,
    reactant_smiles: List[str],
) -> Dict[int, Dict[str, Any]]:
    """Map atom map numbers to element/aromatic/reactant index."""
    if not mapped_smiles or ">>" not in mapped_smiles:
        return {}
    reactants_side = mapped_smiles.split(">>", 1)[0]
    components = [c for c in reactants_side.split(".") if c]

    # Canonical reactant smiles for matching
    canon_reactants: List[str] = []
    for smi in reactant_smiles:
        canon = rdkit_helpers.canonical_smiles(smi) or smi
        canon_reactants.append(canon)
    remaining: Dict[str, List[int]] = {}
    for idx, smi in enumerate(canon_reactants):
        remaining.setdefault(smi, []).append(idx)

    info: Dict[int, Dict[str, Any]] = {}
    for comp in components:
        comp_unmapped = _strip_atom_mapping(comp)
        comp_canon = rdkit_helpers.canonical_smiles(comp_unmapped) or comp_unmapped
        react_idx = None
        if comp_canon in remaining and remaining[comp_canon]:
            react_idx = remaining[comp_canon].pop(0)
        mol = rdkit_helpers.parse_smiles(comp)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if not map_num:
                continue
            info[map_num] = {
                "element": atom.GetSymbol(),
                "aromatic": bool(atom.GetIsAromatic()),
                "reactant_idx": react_idx,
            }
    return info


def _find_halide_for_carbon(
    carbon_map: int,
    broken_bonds: Iterable[Any],
    map_info: Dict[int, Dict[str, Any]],
) -> Optional[str]:
    priority = {"I": 3, "Br": 2, "Cl": 1, "F": 0}
    found: List[str] = []
    for bond in broken_bonds or []:
        if not isinstance(bond, (tuple, list)) or len(bond) < 2:
            continue
        a, b = bond[0], bond[1]
        if a == carbon_map:
            other = b
        elif b == carbon_map:
            other = a
        else:
            continue

        if isinstance(other, int):
            element = map_info.get(other, {}).get("element")
        else:
            element = _parse_element(str(other))
        if element in priority:
            found.append(element)
    if not found:
        return None
    found.sort(key=lambda el: -priority.get(el, 0))
    return found[0]


def _reacted_motifs_by_reactant(
    reactant_bundles: Iterable[Dict[str, Any]],
    reacted_set: Set[str],
) -> List[List[str]]:
    buckets: List[List[str]] = []
    for bundle in reactant_bundles or []:
        motifs = extract_motif_ids(bundle.get("motifs", []))
        buckets.append([m for m in motifs if m in reacted_set])
    return buckets


def _pair_reactants_from_mapping(
    reacted_set: Set[str],
    *,
    analysis: Optional[Dict[str, Any]],
    reactant_smiles: List[str],
    reactant_bundles: List[Dict[str, Any]],
) -> List[str]:
    if not analysis:
        return []
    mapped_smiles = analysis.get("mapped_smiles") or ""
    broken_bonds = analysis.get("broken_bonds") or []
    formed_bonds = analysis.get("formed_bonds") or []
    if not mapped_smiles or not formed_bonds:
        return []

    map_info = _build_map_info(mapped_smiles, reactant_smiles)
    if not map_info:
        return []

    reacted_by_idx = _reacted_motifs_by_reactant(reactant_bundles, reacted_set)
    pairs: List[str] = []

    for bond in formed_bonds:
        if not isinstance(bond, (tuple, list)) or len(bond) < 2:
            continue
        a, b = bond[0], bond[1]
        if not isinstance(a, int) or not isinstance(b, int):
            continue
        a_info = map_info.get(a)
        b_info = map_info.get(b)
        if not a_info or not b_info:
            continue
        # Identify hetero vs carbon.
        if a_info["element"] == "C" and b_info["element"] in {"S", "N", "O"}:
            carbon_map, hetero_map = a, b
            hetero_el = b_info["element"]
        elif b_info["element"] == "C" and a_info["element"] in {"S", "N", "O"}:
            carbon_map, hetero_map = b, a
            hetero_el = a_info["element"]
        else:
            continue

        nuc_idx = map_info.get(hetero_map, {}).get("reactant_idx")
        elec_idx = map_info.get(carbon_map, {}).get("reactant_idx")
        if nuc_idx is None or elec_idx is None:
            continue
        if nuc_idx >= len(reacted_by_idx) or elec_idx >= len(reacted_by_idx):
            continue

        nuc_candidates = reacted_by_idx[nuc_idx]
        elec_candidates = reacted_by_idx[elec_idx]
        if not nuc_candidates or not elec_candidates:
            continue

        if hetero_el == "S":
            nuc = next((m for m in nuc_candidates if m.endswith("-SH")), nuc_candidates[0])
        elif hetero_el == "N":
            nuc = next((m for m in nuc_candidates if m.endswith("-NH2") or m.endswith("-NHR")), nuc_candidates[0])
        elif hetero_el == "O":
            nuc = next((m for m in nuc_candidates if m.endswith("-OH")), nuc_candidates[0])
        else:
            nuc = nuc_candidates[0]

        halide = _find_halide_for_carbon(carbon_map, broken_bonds, map_info)
        if halide:
            elec = next((m for m in elec_candidates if m.endswith(f"-{halide}")), elec_candidates[0])
        else:
            elec = elec_candidates[0]

        if nuc and elec:
            pairs.append(f"{nuc}~{elec}")

    out: List[str] = []
    seen: Set[str] = set()
    for p in pairs:
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def _fallback_pairs_by_halide_priority(reacted: Iterable[str]) -> List[str]:
    """Fallback pairing when mapping is unavailable."""
    reacted_list = [str(r) for r in reacted if r]
    if not reacted_list:
        return []
    nucleophile = next((m for m in reacted_list if m.endswith("-SH")), None)
    if not nucleophile:
        nucleophile = next((m for m in reacted_list if m.endswith("-NH2") or m.endswith("-NHR")), None)
    if not nucleophile:
        nucleophile = next((m for m in reacted_list if m.endswith("-OH")), None)
    if not nucleophile:
        return []
    halide_priority = ["-I", "-Br", "-Cl", "-F"]
    electrophile = None
    for h in halide_priority:
        electrophile = next((m for m in reacted_list if m.endswith(h)), None)
        if electrophile:
            break
    if not electrophile:
        return []
    return [f"{nucleophile}~{electrophile}"]


def _infer_product_broad_tags(bond_key: Optional[str]) -> List[str]:
    tags: List[str] = []
    formed = _extract_bond_section(bond_key, section="form")
    if not formed:
        return tags
    for token in formed:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        if (_is_carbon_label(left) and _is_sulfur_label(right)) or (
            _is_sulfur_label(left) and _is_carbon_label(right)
        ):
            tags.append("Product_C-S")
            if left.endswith("(ar)") or right.endswith("(ar)"):
                tags.append("Product_Aryl_S")
            break
    return sorted(set(tags))


def _select_primary_broad_tag(tags: Iterable[str]) -> str:
    if not tags:
        return "[]"
    priority = {
        "Product_Aryl_S": 2,
        "Product_C-S": 1,
    }
    candidates = [str(t) for t in tags if t]
    if not candidates:
        return "[]"
    candidates.sort(key=lambda t: (-priority.get(t, 0), t))
    return candidates[0]


def _detect_aryl_s_from_product_smiles(product_smiles: Iterable[str]) -> bool:
    if not rdkit_helpers.rdkit_available():
        return False
    # Aromatic carbon bonded to sulfur (aryl thioether/aryl sulfur link)
    # Matches c-S- and c-S(=O)- variants loosely.
    smarts = "[c][S]"
    try:
        from chemtools.util.smarts_cache import compile_smarts
    except Exception:
        return False
    patt = compile_smarts(smarts, validate=False)
    if patt is None:
        return False
    for smi in product_smiles or []:
        mol = rdkit_helpers.parse_smiles(smi)
        if mol is None:
            continue
        try:
            if mol.HasSubstructMatch(patt):
                return True
        except Exception:
            continue
    return False


def _infer_product_broad_tags_with_validation(
    *,
    bond_key: Optional[str],
    product_smiles: Iterable[str],
) -> List[str]:
    """Infer broad product tags with mapping-first + SMARTS fallback."""
    tags = _infer_product_broad_tags(bond_key)
    if "Product_Aryl_S" in tags:
        return tags
    # Fallback: product structure contains aryl-S bond.
    if _detect_aryl_s_from_product_smiles(product_smiles):
        tags.append("Product_C-S")
        tags.append("Product_Aryl_S")
    return sorted(set(tags))


def format_crk_key(
    *,
    bond_key: Optional[str],
    reacted: Iterable[str],
    spectators: Iterable[str],
    product_broad_tags: Iterable[str],
    product_motifs_reactive: Optional[Iterable[str]] = None,
    include_product: bool = False,
) -> str:
    """Build a composite Condition Recommendation Key (CRK-v1)."""
    reactants_text = _format_motif_list(reacted)
    broad_tags = sorted(str(t) for t in product_broad_tags if t)
    reactive_products = sorted(str(t) for t in (product_motifs_reactive or []) if t)

    products_primary = "[]"
    if include_product:
        if reactive_products:
            products_primary = "|".join(reactive_products)
        else:
            products_primary = _select_primary_broad_tag(broad_tags)

    summary = f"|{reactants_text} -> {products_primary}"

    sections: List[str] = [summary]

    if bond_key:
        formed = _extract_bond_section(bond_key, section="form")
        broken = _extract_bond_section(bond_key, section="break")
        if formed:
            sections.append("bond_formed: " + "; ".join(formed))
            if product_motifs_reactive:
                labeled = _label_formed_bonds(formed, product_motifs_reactive)
                if labeled:
                    sections.append("bond_formed_labeled: " + "; ".join(labeled))
        if broken:
            sections.append("bond_broken: " + "; ".join(broken))

    if spectators:
        sections.append("spectators: " + _format_motif_list(spectators))

    return " | ".join(sections)


def _format_motif_list(items: Iterable[str]) -> str:
    values = sorted(str(item) for item in items if item)
    return "|".join(values) if values else "[]"


def _project_formed_motifs_by_taxonomy(
    *,
    reaction_type: Optional[str],
    formed_in_product: Set[str],
    inferred_in_product: Iterable[str],
) -> List[str]:
    allowed = _taxonomy_allowed_product_motifs(reaction_type)
    if not allowed:
        return []
    projected = set(formed_in_product) | {str(m) for m in inferred_in_product if m}
    projected &= allowed
    return sorted(projected)


def _target_product_motifs_from_formed_bonds(
    *,
    formed_bonds: Iterable[str],
    formed_set: Set[str],
) -> Set[str]:
    logic_sets = _load_compound_logic_sets()
    target_ids: Set[str] = set()

    def add_set(name: str) -> None:
        target_ids.update(logic_sets.get(name, set()))

    for token in formed_bonds:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        left_ar = _is_aromatic_label(left)
        right_ar = _is_aromatic_label(right)

        elements = {left_el, right_el}
        if "C" in elements and "N" in elements:
            add_set("aryl_amines")
            add_set("aryl_amides")
            add_set("amides")
        if "C" in elements and "O" in elements:
            add_set("aryl_ethers")
        if "C" in elements and "S" in elements:
            add_set("aryl_thioethers")
        if left_el == "C" and right_el == "C" and left_ar and right_ar:
            target_ids.add("Ar-Ar")
        if left_el == "C" and right_el == "C" and (left_ar ^ right_ar):
            for candidate in ("Ar-Alkyl", "Ar-Alkenyl", "Ar-Alkynyl"):
                if candidate in formed_set:
                    target_ids.add(candidate)

    return target_ids


def _select_reactive_product_motifs(
    product_motifs: Iterable[Dict[str, Any]],
    *,
    bond_key: Optional[str],
    formed_motifs: Iterable[str],
    reacted_motifs: Iterable[str],
    reaction_type: Optional[str] = None,
) -> List[str]:
    """Select reaction-center product motifs with taxonomy constraints first."""
    motif_ids = [
        normalize_motif_id(str(m.get("compound_id") or m.get("id")))
        for m in product_motifs
        if isinstance(m, dict) and (m.get("compound_id") or m.get("id"))
    ]
    motif_set = set(motif_ids)
    formed_set = {normalize_motif_id(str(m)) for m in formed_motifs if m}
    reacted_set = {normalize_motif_id(str(m)) for m in reacted_motifs if m}
    formed_in_product = (motif_set & formed_set) if motif_set else set(formed_set)

    inferred = _infer_product_motifs_from_logic(reacted_set, bond_key)
    inferred_norm = [normalize_motif_id(str(m)) for m in inferred if m]
    inferred_in_product = [m for m in inferred_norm if m in motif_set] if motif_set else inferred_norm
    rescue_without_bond = _rescue_decarboxylative_product_motifs_without_bond_key(
        reacted_motifs=reacted_set,
        product_motif_ids=motif_set or formed_set,
    )
    if rescue_without_bond:
        inferred_in_product = sorted(set(inferred_in_product) | set(rescue_without_bond))

    taxonomy_projected = _project_formed_motifs_by_taxonomy(
        reaction_type=reaction_type,
        formed_in_product=formed_in_product,
        inferred_in_product=inferred_in_product,
    )
    if taxonomy_projected:
        return taxonomy_projected

    if not bond_key:
        if inferred_in_product:
            return inferred_in_product
        if formed_in_product:
            return sorted(formed_in_product)
        if formed_set:
            return sorted(formed_set)
        return []

    formed_bonds = _extract_bond_section(bond_key, section="form")
    if not formed_bonds:
        if inferred_in_product:
            return inferred_in_product
        if formed_in_product:
            return sorted(formed_in_product)
        if formed_set:
            return sorted(formed_set)
        return []

    target_ids = _target_product_motifs_from_formed_bonds(
        formed_bonds=formed_bonds,
        formed_set=formed_set,
    )
    if motif_set:
        reactive = sorted((motif_set & target_ids) & formed_in_product)
    else:
        reactive = sorted(formed_set & target_ids)
    if inferred_in_product:
        if reactive:
            return sorted(set(reactive) | set(inferred_in_product))
        return inferred_in_product
    if reactive:
        return reactive
    if formed_in_product:
        return sorted(formed_in_product)
    if formed_set:
        return sorted(formed_set)
    return []


def _rescue_decarboxylative_product_motifs_without_bond_key(
    *,
    reacted_motifs: Iterable[str],
    product_motif_ids: Iterable[str],
) -> List[str]:
    """
    Rescue decarboxylative C-H functionalization product motifs without bond_key.

    This is used when mapping-derived bond deltas are unavailable/unreliable.
    """
    reacted_set = {normalize_motif_id(str(m)) for m in reacted_motifs if m}
    product_set = {normalize_motif_id(str(m)) for m in product_motif_ids if m}
    if not reacted_set or not product_set:
        return []

    aromatic_ch_sites = {"Ar-H", "HeteroAr-H", "AromN-H"}
    has_aromatic_ch_substrate = any(m in aromatic_ch_sites for m in reacted_set)
    if not has_aromatic_ch_substrate:
        return []

    logic_sets = _load_compound_logic_sets()
    carboxylic_acids = logic_sets.get("carboxylic_acids", set())
    has_carboxylate_partner = any(
        (m in carboxylic_acids) or m.endswith("-CO2H") or m.endswith("-COOH")
        for m in reacted_set
    )
    if not has_carboxylate_partner:
        return []

    candidates: Set[str] = set()
    # Acylation classes
    if "Ar-COR" in product_set:
        candidates.add("Ar-COR")
    if any(m in {"HeteroAr-H", "AromN-H"} for m in reacted_set) and "HeteroAr-COR" in product_set:
        candidates.add("HeteroAr-COR")
    # Alkylation classes
    for motif_id in ("Ar-Alkyl", "Ar-Alkenyl", "Ar-Alkynyl"):
        if motif_id in product_set:
            candidates.add(motif_id)

    return sorted(candidates)


def _infer_product_motifs_from_logic(
    reacted_motifs: Iterable[str],
    bond_key: Optional[str],
) -> List[str]:
    """Infer likely product motifs from reacted motifs + bond formation."""
    if not bond_key:
        return []
    reacted_list = [str(m) for m in reacted_motifs if m]
    if not reacted_list:
        return []
    formed_bonds = _extract_bond_section(bond_key, section="form")
    if not formed_bonds:
        return []
    broken_bonds = _extract_bond_section(bond_key, section="break")

    logic_sets = _load_compound_logic_sets()
    amines = logic_sets.get("aryl_amines", set())
    amides = logic_sets.get("aryl_amides", set())
    ethers = logic_sets.get("aryl_ethers", set())
    thioethers = logic_sets.get("aryl_thioethers", set())
    inference_rules = _load_group_inference_rules()
    suffix_rules = inference_rules.get("suffixes", {}) or {}
    thiol_suffixes = tuple(str(s) for s in (suffix_rules.get("thiol", []) or []) if s)
    alcohol_suffixes = tuple(str(s) for s in (suffix_rules.get("alcohol", []) or []) if s)
    amide_suffixes = tuple(str(s) for s in (suffix_rules.get("amide", []) or []) if s)
    amine_primary_suffixes = tuple(str(s) for s in (suffix_rules.get("amine_primary", []) or []) if s)
    amine_secondary_suffixes = tuple(str(s) for s in (suffix_rules.get("amine_secondary", []) or []) if s)
    amine_tertiary_suffixes = tuple(str(s) for s in (suffix_rules.get("amine_tertiary", []) or []) if s)
    aromatic_nitrogen = {
        str(s) for s in (inference_rules.get("aromatic_nitrogen", []) or []) if s
    }
    alkyl_prefixes = tuple(str(s) for s in (inference_rules.get("alkyl_prefixes", []) or []) if s)
    alkyl_leaving = tuple(str(s) for s in (inference_rules.get("alkyl_leaving", []) or []) if s)
    alkenyl_leaving = set(str(s) for s in (inference_rules.get("alkenyl_leaving", []) or []) if s)
    alkynyl_leaving = set(str(s) for s in (inference_rules.get("alkynyl_leaving", []) or []) if s)

    def has_suffix(suffixes: Tuple[str, ...]) -> bool:
        return any(m.endswith(sfx) for m in reacted_list for sfx in suffixes)

    def has_any(values: Set[str]) -> bool:
        if not values:
            return False
        for m in reacted_list:
            if m in values:
                return True
        return False

    inferred: Set[str] = set()
    for token in formed_bonds:
        if "-" not in token:
            continue
        left, right = [t.strip() for t in token.split("-", 1)]
        left_el = _parse_element(left)
        right_el = _parse_element(right)
        left_ar = _is_aromatic_label(left)
        right_ar = _is_aromatic_label(right)
        elements = {left_el, right_el}
        is_aryl_carbon = (left_el == "C" and left_ar) or (right_el == "C" and right_ar)
        is_c_c = left_el == "C" and right_el == "C"
        is_aryl_sp3_c = is_c_c and (left_ar ^ right_ar)

        if "C" in elements and "S" in elements and is_aryl_carbon:
            if has_any(logic_sets.get("thiols_sh", set())) or has_suffix(thiol_suffixes):
                if "Ar-SR" in thioethers:
                    inferred.add("Ar-SR")

        if "C" in elements and "O" in elements and is_aryl_carbon:
            if has_any(logic_sets.get("alcohols_oh", set())) or has_suffix(alcohol_suffixes):
                if "Ar-OR" in ethers:
                    inferred.add("Ar-OR")

        if "C" in elements and "N" in elements and is_aryl_carbon:
            n_is_aromatic = (left_el == "N" and left_ar) or (right_el == "N" and right_ar)
            if n_is_aromatic:
                if "Ar-AromN" in amines:
                    inferred.add("Ar-AromN")
                continue

            if has_suffix(amide_suffixes):
                if "Ar-NRCOR" in amides:
                    inferred.add("Ar-NRCOR")
                continue
            if has_suffix(amine_primary_suffixes):
                if "Ar-NH2" in amines:
                    inferred.add("Ar-NH2")
                continue
            if has_suffix(amine_secondary_suffixes):
                if "Ar-NHR" in amines:
                    inferred.add("Ar-NHR")
                continue
            if has_suffix(amine_tertiary_suffixes):
                if "Ar-NR2" in amines:
                    inferred.add("Ar-NR2")
                continue

            if aromatic_nitrogen and any(m in aromatic_nitrogen for m in reacted_list) and "Ar-AromN" in amines:
                inferred.add("Ar-AromN")

        if is_aryl_sp3_c:
            has_alkyl_electrophile = False
            if alkyl_prefixes and alkyl_leaving:
                has_alkyl_electrophile = any(
                    m.startswith(alkyl_prefixes) and m.endswith(alkyl_leaving) for m in reacted_list
                )
            if has_alkyl_electrophile:
                inferred.add("Ar-Alkyl")
            if alkenyl_leaving and (
                any(
                    m in alkenyl_leaving or (alkyl_leaving and m.startswith("Alkenyl-") and m.endswith(alkyl_leaving))
                    for m in reacted_list
                )
            ):
                inferred.add("Ar-Alkenyl")
            if alkynyl_leaving and (
                any(
                    m in alkynyl_leaving or (alkyl_leaving and m.startswith("Alkynyl-") and m.endswith(alkyl_leaving))
                    for m in reacted_list
                )
            ):
                inferred.add("Ar-Alkynyl")

    # Decarboxylative aromatic C-H functionalization (Minisci-like):
    # infer product class from bond topology + reactant motifs, even when the
    # same product motif family is present on both sides and delta counts are zero.
    has_aryl_c_c_formation = any(
        "-" in token
        and _parse_element(token.split("-", 1)[0].strip()) == "C"
        and _parse_element(token.split("-", 1)[1].strip()) == "C"
        and (
            _is_aromatic_label(token.split("-", 1)[0].strip())
            or _is_aromatic_label(token.split("-", 1)[1].strip())
        )
        for token in formed_bonds
    )
    has_c_c_cleavage = any(
        "-" in token
        and _parse_element(token.split("-", 1)[0].strip()) == "C"
        and _parse_element(token.split("-", 1)[1].strip()) == "C"
        for token in broken_bonds
    )
    aromatic_ch_sites = {"Ar-H", "HeteroAr-H", "AromN-H"}
    has_aromatic_ch_substrate = any(m in aromatic_ch_sites for m in reacted_list)
    carboxylic_acids = logic_sets.get("carboxylic_acids", set())
    has_carboxylate_partner = any(
        (m in carboxylic_acids) or m.endswith("-CO2H") or m.endswith("-COOH")
        for m in reacted_list
    )
    has_acyl_ketoacid_partner = any(
        m == "Acyl-CO2H" or (m.startswith("Acyl-") and m.endswith("-CO2H"))
        for m in reacted_list
    )
    if has_aryl_c_c_formation and has_c_c_cleavage and has_aromatic_ch_substrate and has_carboxylate_partner:
        if has_acyl_ketoacid_partner:
            # Keep Ar-COR for broad taxonomy compatibility and add heteroaryl
            # variant when the reacting C-H site is heteroaromatic.
            inferred.add("Ar-COR")
            if any(m in {"HeteroAr-H", "AromN-H"} for m in reacted_list):
                inferred.add("HeteroAr-COR")
        else:
            inferred.add("Ar-Alkyl")

    return sorted(inferred)


def _scaffold_spectators_from_bundles(
    reactant_bundles: Iterable[Dict[str, Any]],
    product_bundles: Iterable[Dict[str, Any]],
) -> Set[str]:
    """Identify scaffold motifs present on both sides to treat as spectators in CRK."""
    try:
        from .aggregation import load_scaffold_motif_ids
    except Exception:
        return set()
    scaffold_ids = load_scaffold_motif_ids()
    if not scaffold_ids:
        return set()
    reactant_ids: Set[str] = set()
    product_ids: Set[str] = set()
    for bundle in reactant_bundles or []:
        reactant_ids.update(
            extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", []))
        )
    for bundle in product_bundles or []:
        product_ids.update(
            extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", []))
        )
    return set(m for m in reactant_ids & product_ids if m in scaffold_ids)


def classify_agent_roles(agents: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    """Classify reagents/solvents by role using reagent taxonomy."""
    from functools import lru_cache
    from chemtools.reagent.reagent_v2 import ReagentTaxonomyV2
    
    @lru_cache(maxsize=1)
    def load_taxonomy() -> Optional[ReagentTaxonomyV2]:
        try:
            return ReagentTaxonomyV2.from_path()
        except Exception:
            return None
    
    def get_agent_smiles(agent: Dict[str, Any]) -> str:
        """Extract SMILES from agent dict."""
        for key in ("smiles", "smiles_norm", "largest_smiles"):
            value = agent.get(key)
            if value:
                return str(value)
        return ""
    
    taxonomy = load_taxonomy()
    if not taxonomy:
        return {"agent_count": 0, "role_counts": {}, "family_counts": {}, "role_flags": {}, "flags": {}}
    
    role_flags = [
        "metal_catalyst", "organo_catalyst", "enzyme", "ligand", "base",
        "acid", "solvent", "additive", "oxidant", "reductant",
        "condensation_agent", "other_reagent"
    ]
    
    entries: List[Dict[str, Any]] = []
    role_counts: Dict[str, int] = {}
    family_counts: Dict[str, int] = {}
    flags: Dict[str, bool] = {role: False for role in role_flags}
    
    for agent in agents or []:
        smiles = get_agent_smiles(agent)
        if not smiles:
            continue
        
        reagent = taxonomy.lookup_reagent(smiles)
        if not reagent:
            continue
        
        role = reagent.get("role") or "other_reagent"
        family = reagent.get("family") or "Unknown"
        
        role_counts[role] = role_counts.get(role, 0) + 1
        family_counts[family] = family_counts.get(family, 0) + 1
        
        if role in flags:
            flags[role] = True
        
        entries.append({
            "smiles": smiles,
            "role": role,
            "family": family,
            "name": reagent.get("name"),
        })
    
    return {
        "entries": entries,
        "agent_count": len(entries),
        "role_counts": role_counts,
        "family_counts": family_counts,
        "role_flags": {k: v for k, v in flags.items() if v},
        "flags": flags,
    }


def featurize_reaction(
    reaction_smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Return a reaction feature bundle (core or extended format).
    
    Args:
        reaction_smiles: Reaction SMILES with >> separator
        registry_paths: Custom taxonomy paths
        options: Featurization options
            - detailed (bool): If True, return extended format with all analysis.
                             If False (default), return simplified core format.
        
    Returns:
        Core bundle (7 fields) or extended bundle (8 fields with extended section)
    """
    options = options or {}
    use_detailed = to_bool(options.get("detailed"), default=False)
    
    include_roles = to_bool(options.get("include_roles"), default=True)
    include_agent_roles = to_bool(options.get("include_agent_roles"), default=True)
    include_reaction_type = to_bool(options.get("include_reaction_type"), default=False)
    skip_bond_analysis = to_bool(options.get("skip_bond_analysis"), default=False)
    include_product_in_crk = to_bool(options.get("include_product_in_crk"), default=True)
    reactant_coverage_guard = to_bool(options.get("reactant_coverage_guard"), default=True)
    
    # General coupling confirmation (supports 9+ coupling reaction types)
    # Backward compatibility: map old Suzuki-specific parameter to general one
    confirm_coupling = options.get("confirm_coupling_products")
    if confirm_coupling is None:
        # Check for legacy Suzuki-specific parameter
        confirm_coupling = options.get("confirm_suzuki_products")
    confirm_coupling_products = to_bool(confirm_coupling, default=True)

    # Start with Unknown, then determine from CRK_raw after bond/motif extraction.
    reaction_type = {"reaction_type": "Unknown", "confidence": 0.0, "slot_evidence": {}}
    detection_payload = {}
    quality_warnings: List[str] = []

    # Normalize reaction SMILES
    normalized = normalize_reaction(reaction_smiles)
    reaction_record = ReactionRecord.from_payload(normalized)
    reactant_smiles = reaction_record.reactant_smiles

    if include_reaction_type:
        primary_detection = _run_primary_reaction_type_detection(reaction_smiles)
        detection_payload["primary_detection"] = {
            "status": primary_detection.get("status"),
        }
        primary_error = primary_detection.get("error")
        if primary_error:
            detection_payload["primary_detection"]["error"] = primary_error
        primary_matches = primary_detection.get("matches") or []
        if isinstance(primary_matches, list) and primary_matches:
            detection_payload["matches"] = primary_matches
        primary_summary = primary_detection.get("summary")
        if isinstance(primary_summary, dict):
            rt_candidate = _resolve_reaction_type_id(primary_summary.get("reaction_type"))
            rt_conf = _to_float_or_default(primary_summary.get("confidence"), 0.0)
            if rt_candidate:
                reaction_type = _set_reaction_type_payload(
                    reaction_type,
                    rt_candidate,
                    rt_conf,
                )

    # Classify agents/reagents
    agent_roles = None
    if include_agent_roles:
        agent_roles = classify_agent_roles(reaction_record.agent_payloads)

    # Featurize reactants
    reactant_bundles = [
        build_molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        for smiles in reactant_smiles
    ]
    _ensure_reactant_coverage(reactant_bundles, enabled=reactant_coverage_guard)

    # Featurize products
    product_smiles = reaction_record.product_smiles
    product_bundles: List[Dict[str, Any]] = []
    product_motif_ids: List[str] = []
    product_motifs_full: List[Dict[str, Any]] = []  # Full motif dicts with fingerprints
    for smiles in product_smiles:
        try:
            bundle = build_molecule_bundle(smiles, registry_paths=registry_paths, options=options)
        except Exception:
            continue
        product_bundles.append(bundle)
        product_motif_ids.extend(extract_motif_ids(bundle.get("motifs", []), bundle.get("context_motifs", [])))
        # Collect full motif dicts for fingerprint-aware comparison
        product_motifs_full.extend(bundle.get("motifs", []))
        product_motifs_full.extend(bundle.get("context_motifs", []))

    # Aggregate features without pattern-based filtering
    # Pattern filtering is skipped to avoid removing reacted motifs based on incorrect initial detection
    # The validation step will correct the detection using unfiltered reacted motifs
    aggregates = aggregate_reaction_features(
        reactant_bundles,
        product_motif_ids=product_motif_ids,
        product_motifs=product_motifs_full,
        reaction_type=None,  # Disable pattern-based filtering
    )
    
    # Classify reactant roles
    roles_summary = None
    if include_roles:
        try:
            roles = classify_reactants_with_context(reaction_smiles)
            roles_summary = get_reactant_summary(roles)
        except Exception:
            roles_summary = None

    intramolecular = infer_intramolecular(reactant_smiles, product_smiles, roles_summary)

    rt_id = None
    llm_assist_meta: Optional[Dict[str, Any]] = None
    reaction_events: Dict[str, Any] = {}

    # Generate CRK-v1 reaction key (single source of truth)
    reaction_key = None
    product_broad_tags: List[str] = []
    product_motifs_reactive: List[str] = []
    if product_bundles and reactant_bundles:
        reacted_full = set(aggregates.get("reacted_motifs", []))
        spectators = set(aggregates.get("spectator_motifs", []))
        scaffold_spectators = _scaffold_spectators_from_bundles(reactant_bundles, product_bundles)
        spectators_for_crk = spectators | scaffold_spectators
        if scaffold_spectators:
            group_list = list(aggregates.get("spectator_groups_combined") or [])
            seen_groups = {str(g).strip() for g in group_list if str(g).strip()}
            for scaffold_id in sorted(scaffold_spectators):
                sid = str(scaffold_id).strip()
                if not sid or sid in seen_groups:
                    continue
                seen_groups.add(sid)
                group_list.append(sid)
            aggregates["spectator_groups_combined"] = group_list
            aggregates["spectator_groups_ranked"] = rank_spectator_groups(group_list)

        bond_analysis = None
        fallback_bond_analysis = None
        bond_key = None
        fallback_bond_key = None
        if not skip_bond_analysis:
            preferred_bond_analysis, mapping_unreliable, mapping_meta = _get_bond_change_analysis_with_quality(
                reaction_smiles
            )
            if preferred_bond_analysis:
                if mapping_unreliable:
                    fallback_bond_analysis = preferred_bond_analysis
                    agreement_payload = (mapping_meta or {}).get("agreement") or {}
                    detection_payload["mapping_warning"] = {
                        "reason": "low_confidence_mapping_disagreement",
                        "fallback_used_for_product_projection": True,
                        "combined_confidence": _to_float_or_default(
                            (mapping_meta or {}).get("combined_confidence"),
                            0.0,
                        ),
                        "agreement_rxnmapper_vs_mcs": agreement_payload.get("rxnmapper_vs_mcs"),
                        "validation": (mapping_meta or {}).get("validation"),
                    }
                else:
                    bond_analysis = preferred_bond_analysis
            if bond_analysis:
                bond_key = format_bond_change_key(reaction_smiles, analysis=bond_analysis)
            if fallback_bond_analysis:
                fallback_bond_key = format_bond_change_key(
                    reaction_smiles,
                    analysis=fallback_bond_analysis,
                )
        bond_key = _sanitize_bond_key(
            bond_key,
            reacted_full,
            group_element_map=_load_group_element_map(),
        )
        fallback_bond_key = _sanitize_bond_key(
            fallback_bond_key,
            reacted_full,
            group_element_map=_load_group_element_map(),
        )
        projection_bond_key = bond_key or fallback_bond_key

        # Determine reaction type from CRK_raw (before product projection).
        reacted_for_detection = _filter_reactants_for_crk(
            reacted_full,
            bond_key,
            spectators=spectators_for_crk,
        )
        if not reacted_for_detection and reacted_full:
            reacted_for_detection = sorted(
                normalize_motif_id(str(m))
                for m in reacted_full
                if m
            )
        formed_raw = {
            normalize_motif_id(str(m))
            for m in (aggregates.get("formed_motifs") or [])
            if m
        }
        formed_inferred = {
            normalize_motif_id(str(m))
            for m in _infer_product_motifs_from_logic(
                reacted_for_detection,
                bond_key,
            )
            if m
        }
        if not formed_inferred:
            rescue_formed = _rescue_decarboxylative_product_motifs_without_bond_key(
                reacted_motifs=reacted_for_detection,
                product_motif_ids=(
                    normalize_motif_id(str(m.get("compound_id") or m.get("id")))
                    for m in product_motifs_full
                    if isinstance(m, dict) and (m.get("compound_id") or m.get("id"))
                ),
            )
            formed_inferred = {normalize_motif_id(str(m)) for m in rescue_formed if m}
        formed_all = sorted(formed_raw | formed_inferred)
        spectators_for_detection = sorted(
            normalize_motif_id(str(m))
            for m in (set(spectators_for_crk) - set(reacted_for_detection))
            if m
        )
        reaction_key_raw = format_crk_key(
            bond_key=bond_key,
            reacted=reacted_for_detection,
            spectators=spectators_for_detection,
            product_broad_tags=[],
            product_motifs_reactive=formed_all,
            include_product=True,
        )
        from .detection_validation import validate_detection_with_crk_key

        validated = validate_detection_with_crk_key(
            initial_detection=(
                reaction_type.get("reaction_type", "Unknown")
                if isinstance(reaction_type, dict)
                else str(reaction_type)
            ),
            initial_confidence=(
                reaction_type.get("confidence", 0.0)
                if isinstance(reaction_type, dict)
                else 0.0
            ),
            reaction_key=reaction_key_raw,
        )
        if validated.get("reaction_type"):
            if isinstance(reaction_type, dict):
                reaction_type["reaction_type"] = validated["reaction_type"]
                reaction_type["name"] = validated["reaction_type"]
                reaction_type["confidence"] = validated["confidence"]
            else:
                reaction_type = validated["reaction_type"]
        detection_payload["validation"] = {
            "original_detection": validated.get("corrected_from"),
            "validated_detection": validated.get("reaction_type"),
            "validation_method": validated.get("validation_method"),
            "validation_reason": validated.get("reason"),
            "validation_confidence": validated.get("confidence"),
            "reaction_key_raw": reaction_key_raw,
        }
        if validated.get("evidence") is not None:
            detection_payload["evidence"] = validated.get("evidence")

        if isinstance(reaction_type, dict):
            rt_id = reaction_type.get("reaction_type")
        elif reaction_type is not None:
            rt_id = str(reaction_type)
        if rt_id == "Unknown":
            rt_id = None

        product_broad_tags = _infer_product_broad_tags_with_validation(
            bond_key=projection_bond_key,
            product_smiles=product_smiles,
        )
        product_motifs_reactive = _select_reactive_product_motifs(
            product_motifs_full,
            bond_key=projection_bond_key,
            formed_motifs=formed_all,
            reacted_motifs=reacted_for_detection,
            reaction_type=rt_id,
        )
        formed_center = sorted(set(product_motifs_reactive))
        formed_center_set = set(formed_center)
        formed_context = sorted(m for m in formed_all if m not in formed_center_set)
        aggregates["formed_motifs_all"] = formed_all
        aggregates["formed_motifs_center"] = formed_center
        aggregates["formed_motifs_context"] = formed_context
        reacted_for_crk = list(reacted_for_detection)
        spectators_for_key = sorted(set(spectators_for_crk) - set(reacted_for_crk))
        reaction_key = format_crk_key(
            bond_key=bond_key,
            reacted=reacted_for_crk,
            spectators=spectators_for_key,
            product_broad_tags=product_broad_tags,
            product_motifs_reactive=product_motifs_reactive,
            include_product=include_product_in_crk,
        )

        reaction_events = summarize_reaction_events(
            reaction_smiles=reaction_smiles,
            bond_key=bond_key,
            fallback_bond_key=fallback_bond_key,
            reacted_motifs=aggregates.get("reacted_motifs") or [],
            formed_motifs=aggregates.get("formed_motifs_all")
            or aggregates.get("formed_motifs")
            or [],
            mapping_warning=detection_payload.get("mapping_warning")
            if isinstance(detection_payload, dict)
            else None,
        )
        if reaction_events:
            detection_payload["reaction_key_quality"] = (
                reaction_events.get("reaction_key_quality") or {}
            )
            quality_payload = reaction_events.get("reaction_key_quality") or {}
            quality_level = str(quality_payload.get("level") or "").strip().lower()
            quality_score = _to_float_or_default(quality_payload.get("score_0_1"), 1.0)
            if quality_level == "low" or quality_score < 0.45:
                quality_warnings.append("low_reaction_key_quality")
            anomalies = reaction_events.get("anomalies") or []
            if isinstance(anomalies, list):
                critical_anomalies = {
                    "mapping_unreliable_fallback_used",
                    "amidation_without_explicit_activation_marker",
                }
                for anomaly in anomalies:
                    anomaly_text = str(anomaly).strip()
                    if not anomaly_text:
                        continue
                    if quality_level == "high" and anomaly_text not in critical_anomalies:
                        continue
                    quality_warnings.append(f"reaction_key_anomaly:{anomaly_text}")

    if rt_id is None:
        if isinstance(reaction_type, dict):
            rt_id = reaction_type.get("reaction_type")
        elif reaction_type is not None:
            rt_id = str(reaction_type)
        if rt_id == "Unknown":
            rt_id = None

    llm_assist = _normalize_llm_assist_options(options)
    if llm_assist.get("enabled"):
        llm_assist_meta = {
            "enabled": True,
            "used": False,
            "status": "skipped",
            "decision": "none",
        }
        if not llm_assist.get("provider") or not llm_assist.get("model"):
            llm_assist_meta["status"] = "config_error"
            llm_assist_meta["decision"] = "missing_provider_or_model"
        else:
            uncertain, uncertainty_reasons = _is_reaction_uncertain_for_llm_assist(
                reaction_type=reaction_type,
                detection_payload=detection_payload,
                reaction_key=reaction_key,
                confidence_threshold=_to_float_or_default(
                    llm_assist.get("confidence_threshold"),
                    _LLM_ASSIST_DEFAULT_CONFIDENCE_THRESHOLD,
                ),
            )
            llm_assist_meta["uncertainty_reasons"] = uncertainty_reasons
            only_on_uncertain = bool(llm_assist.get("only_on_uncertain", True))
            if only_on_uncertain and not uncertain:
                llm_assist_meta["status"] = "skipped_not_uncertain"
                llm_assist_meta["decision"] = "deterministic_kept"
            else:
                llm_assist_meta["used"] = True
                current_rt_id, current_conf = _get_reaction_type_id_and_confidence(reaction_type)
                reaction_key_raw = (
                    (detection_payload.get("validation") or {}).get("reaction_key_raw")
                    if isinstance(detection_payload.get("validation"), dict)
                    else None
                )
                review_context = {
                    "reaction_smiles": reaction_smiles,
                    "normalized": normalized,
                    "deterministic_reaction_type": current_rt_id or "Unknown",
                    "deterministic_confidence": current_conf,
                    "reaction_key_raw": reaction_key_raw or "",
                    "reaction_key": reaction_key or "",
                    "mapping_warning": detection_payload.get("mapping_warning"),
                    "reacted_motifs": aggregates.get("reacted_motifs") or [],
                    "formed_motifs": aggregates.get("formed_motifs_all")
                    or aggregates.get("formed_motifs")
                    or [],
                    "spectator_motifs": aggregates.get("spectator_motifs") or [],
                    "product_broad_tags": product_broad_tags or [],
                    "product_motifs_reactive": product_motifs_reactive or [],
                }
                review = _run_llm_reaction_assist(review_context, llm_assist)
                llm_assist_meta["status"] = str(review.get("status") or "unknown")
                llm_assist_meta["provider"] = review.get("provider") or llm_assist.get(
                    "provider"
                )
                llm_assist_meta["model"] = review.get("model") or llm_assist.get("model")
                if review.get("total_tokens") is not None:
                    llm_assist_meta["total_tokens"] = review.get("total_tokens")
                if review.get("latency_ms") is not None:
                    llm_assist_meta["latency_ms"] = review.get("latency_ms")
                if review.get("error"):
                    llm_assist_meta["error"] = review.get("error")

                analysis = review.get("analysis")
                if review.get("status") == "ok" and isinstance(analysis, dict):
                    suggested_label = str(
                        analysis.get("suggested_reaction_type") or ""
                    ).strip()
                    llm_assist_meta["suggested_reaction_type"] = suggested_label or "Unknown"
                    llm_assist_meta["suggested_confidence"] = _to_float_or_default(
                        analysis.get("confidence"),
                        0.0,
                    )
                    llm_assist_meta["requires_human_review"] = bool(
                        analysis.get("requires_human_review", False)
                    )
                    if analysis.get("uncertainty_flags"):
                        llm_assist_meta["uncertainty_flags"] = list(
                            analysis.get("uncertainty_flags") or []
                        )

                    if not suggested_label or suggested_label.lower() == "unknown":
                        llm_assist_meta["decision"] = "no_valid_suggestion"
                    else:
                        suggested_rt_id = _resolve_reaction_type_id(suggested_label)
                        if not suggested_rt_id:
                            llm_assist_meta["decision"] = "invalid_taxonomy_label"
                        elif suggested_rt_id == current_rt_id:
                            llm_assist_meta["decision"] = "no_change"
                        else:
                            validation_ok = True
                            validation_reason = "validation_bypassed"
                            if bool(llm_assist.get("require_crk_validation", True)):
                                validation_ok, validation_reason = (
                                    _validate_llm_suggested_type_with_crk(
                                        suggested_rt_id,
                                        reaction_key_raw or reaction_key,
                                    )
                                )
                            if validation_ok:
                                llm_conf = _to_float_or_default(
                                    analysis.get("confidence"),
                                    current_conf,
                                )
                                merged_conf = max(current_conf, llm_conf)
                                reaction_type = _set_reaction_type_payload(
                                    reaction_type,
                                    suggested_rt_id,
                                    merged_conf,
                                )
                                rt_id = suggested_rt_id
                                llm_assist_meta["decision"] = "applied"
                                llm_assist_meta["validation"] = validation_reason
                            else:
                                llm_assist_meta["decision"] = "rejected_validation_mismatch"
                                llm_assist_meta["validation"] = validation_reason
                elif llm_assist_meta.get("decision") == "none":
                    llm_assist_meta["decision"] = "llm_review_failed"
        detection_payload["llm_assist"] = {
            "decision": llm_assist_meta.get("decision"),
            "status": llm_assist_meta.get("status"),
            "used": llm_assist_meta.get("used", False),
        }

    reaction = {
        "reaction_smiles": reaction_smiles,
        "normalized": normalized,
        "detection": detection_payload,
        "reaction_type": reaction_type,
        "reactants": reactant_bundles,
        "products": product_bundles,
        "aggregates": aggregates,
        "reaction_key": reaction_key,
        "product_broad_tags": product_broad_tags,
        "product_motifs_reactive": product_motifs_reactive,
        "roles": roles_summary,
        "agent_roles": agent_roles,
        "intramolecular": intramolecular,
    }
    if reaction_events:
        reaction["reaction_events"] = reaction_events

    # Add feasibility analysis for specific reaction types
    if rt_id == "snar_cn" or rt_id == "c_n_cross_coupling":
        reaction["feasibility"] = analyze_snar_feasibility(reaction)

    # Return simplified format (core or extended)
    if use_detailed:
        result = build_extended_reaction(reaction)
    else:
        result = build_core_reaction(reaction)
    
    # Add metadata
    result["kind"] = "reaction"
    result["schema_version"] = "v2"
    
    meta = {
        "rdkit_available": rdkit_helpers.rdkit_available(),
    }
    meta_errors: List[str] = []
    if detection_payload.get("error"):
        meta_errors.append(str(detection_payload["error"]))
    if quality_warnings:
        seen_errors: Set[str] = set()
        for warning in quality_warnings:
            text = str(warning).strip()
            if not text or text in seen_errors:
                continue
            seen_errors.add(text)
            meta_errors.append(text)
    if meta_errors:
        meta["errors"] = meta_errors
    if llm_assist_meta is not None:
        meta["llm_assist"] = llm_assist_meta
    if meta.get("errors") or not meta.get("rdkit_available", True):
        result["meta"] = meta
    elif llm_assist_meta is not None:
        result["meta"] = meta
    
    return result
