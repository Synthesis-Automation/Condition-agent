"""Featurizer-based query analysis for recommendation workflows."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .models import QueryAnalysis, RecommendationRequest
from .utils import pick_electrophile_nucleophile


def _normalize_smiles_list(items: List[Dict[str, Any]]) -> List[str]:
    out: List[str] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        smi = item.get("smiles_norm") or item.get("largest_smiles") or item.get("input")
        if smi:
            out.append(str(smi))
    return out


def _resolve_reaction_type(label: Optional[str]) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    if not label:
        return None, None, None, None
    try:
        from chemtools.taxonomy import reaction_catalog
        rid = reaction_catalog.resolve_reaction_type(str(label))
        if not rid:
            return str(label), None, None, None
        rec = reaction_catalog.get_reaction_type(rid)
        if rec is None:
            return rid, rid, None, None
        return rid, getattr(rec, "id", rid), getattr(rec, "name", None), getattr(rec, "category", None)
    except Exception:
        text = str(label).strip() or None
        return text, None, None, None


def _featurize_query_reaction(reaction_smiles: str) -> Dict[str, Any]:
    from chemtools.featurizers.reaction_path import analyze_reaction_featurization

    return analyze_reaction_featurization(
        reaction_smiles,
        base_options={
            "confirm_coupling_products": True,
            "skip_bond_analysis": True,
        },
        cleanup=True,
    ).bundle


def analyze_recommendation_query(request: RecommendationRequest | str) -> QueryAnalysis:
    """Analyze a recommendation query without loading HTE datasets."""
    if isinstance(request, RecommendationRequest):
        req = request
    else:
        req = RecommendationRequest(reaction_smiles=str(request or ""))

    warnings: List[str] = []
    reaction_input = str(req.reaction_smiles or "").strip()
    analysis = QueryAnalysis(
        reaction_smiles_input=reaction_input,
        requested_reaction_type_filter=req.reaction_type_filter,
    )

    if not reaction_input:
        analysis.warnings = ("empty_reaction_smiles",)
        return analysis

    try:
        from chemtools.core.smiles import normalize_reaction

        payload = normalize_reaction(reaction_input)
    except Exception as exc:
        analysis.warnings = (f"normalize_reaction_failed: {exc}",)
        return analysis

    reactants = _normalize_smiles_list(list(payload.get("reactants") or []))
    agents = _normalize_smiles_list(list(payload.get("agents") or []))
    products = _normalize_smiles_list(list(payload.get("products") or []))
    normalized = str(payload.get("normalized") or reaction_input)

    reactant_a = ""
    reactant_b: Optional[str] = None
    if reactants:
        if len(reactants) >= 2:
            try:
                reactant_a, reactant_b = pick_electrophile_nucleophile(reactants)
            except Exception:
                reactant_a = reactants[0]
                reactant_b = ".".join(reactants[1:]) or None
        else:
            reactant_a = reactants[0]

    analysis.reaction_smiles_normalized = normalized
    analysis.reactants = tuple(reactants)
    analysis.agents = tuple(agents)
    analysis.products = tuple(products)
    analysis.reactant_a_smiles = reactant_a
    analysis.reactant_b_smiles = reactant_b
    analysis.product_smiles = products[0] if products else None

    # Resolve user-provided filter through taxonomy when available.
    if req.reaction_type_filter:
        resolved_filter, _, _, _ = _resolve_reaction_type(req.reaction_type_filter)
        analysis.requested_reaction_type_filter_canonical = resolved_filter

    try:
        feat_payload = _featurize_query_reaction(normalized)
        reaction = feat_payload.get("reaction") if isinstance(feat_payload, dict) and isinstance(feat_payload.get("reaction"), dict) else feat_payload
        if not isinstance(reaction, dict):
            raise ValueError("featurize_reaction returned non-dict payload")

        aggregates = reaction.get("aggregates") if isinstance(reaction.get("aggregates"), dict) else {}
        analysis.reaction_key = str(reaction.get("reaction_key") or "").strip() or None
        analysis.reacted_motifs = tuple(str(x) for x in (aggregates.get("reacted_motifs") or []) if str(x).strip())
        analysis.formed_motifs = tuple(str(x) for x in (aggregates.get("formed_motifs") or []) if str(x).strip())
        analysis.spectator_motifs = tuple(str(x) for x in (aggregates.get("spectator_motifs") or []) if str(x).strip())
        spectator_groups = aggregates.get("spectator_groups_ranked") or aggregates.get("spectator_groups_combined") or []
        analysis.spectator_groups = tuple(str(x) for x in spectator_groups if str(x).strip())

        rxn_type_data = reaction.get("reaction_type", {})
        detected_type: Optional[str]
        conf = 0.0
        if isinstance(rxn_type_data, dict):
            detected_type = str(rxn_type_data.get("reaction_type") or rxn_type_data.get("name") or "").strip() or None
            try:
                conf = float(rxn_type_data.get("confidence", 0.0) or 0.0)
            except Exception:
                conf = 0.0
        else:
            detected_type = str(rxn_type_data or "").strip() or None
            try:
                conf = float(reaction.get("confidence", 0.0) or 0.0)
            except Exception:
                conf = 0.0

        resolved, rid, rname, rcat = _resolve_reaction_type(detected_type)
        analysis.detected_reaction_type = resolved or detected_type
        analysis.detected_reaction_type_id = rid
        analysis.detected_reaction_type_name = rname
        analysis.detected_reaction_type_category = rcat
        analysis.reaction_type_confidence = conf

        analysis.raw_feature_summary = {
            "has_reaction_key": bool(analysis.reaction_key),
            "reacted_motif_count": len(analysis.reacted_motifs),
            "formed_motif_count": len(analysis.formed_motifs),
            "spectator_group_count": len(analysis.spectator_groups),
        }
    except Exception as exc:
        warnings.append(f"featurize_reaction_failed: {exc}")

    analysis.warnings = tuple(warnings)
    return analysis

