"""Build relation-preserving condition precedents for web responses."""

from __future__ import annotations

from typing import Any, Dict, Mapping


def _fallback_contexts(recommendation: Mapping[str, Any]) -> list[Dict[str, Any]]:
    """Retain reaction images for older recommendation payloads.

    References are deliberately not guessed from the legacy independent arrays.
    Those arrays are aggregated differently and do not preserve row-level
    relationships.
    """

    reaction_ids = list(recommendation.get("precedent_reaction_ids") or ())
    reaction_smiles = list(
        recommendation.get("precedent_reaction_smiles") or ()
    )
    return [
        {
            "reaction_id": (
                str(reaction_ids[index]) if index < len(reaction_ids) else ""
            ),
            "observation_id": "",
            "reaction_smiles": str(smiles),
            "reference_id": "",
        }
        for index, smiles in enumerate(reaction_smiles)
        if str(smiles).strip()
    ]


def attach_condition_precedents(
    payload: Dict[str, Any],
    reference_catalog: Mapping[str, Mapping[str, Any]],
    experimental_catalog: Mapping[str, Mapping[str, Any]],
) -> Dict[str, Any]:
    """Attach reaction, citation, and procedure as one nested evidence unit."""

    for recommendation in payload.get("recommendations") or ():
        raw_contexts = recommendation.get("precedent_reaction_contexts") or ()
        contexts = [
            dict(context)
            for context in raw_contexts
            if isinstance(context, Mapping)
        ]
        if not contexts:
            contexts = _fallback_contexts(recommendation)

        precedents = []
        for context in contexts:
            reaction_id = str(context.get("reaction_id") or "")
            observation_id = str(context.get("observation_id") or "")
            reference_id = str(context.get("reference_id") or "")
            reference = reference_catalog.get(reference_id)
            experimental = None
            if observation_id:
                experimental = experimental_catalog.get(
                    f"observation:{observation_id}"
                )
            if experimental is None and reaction_id:
                experimental = experimental_catalog.get(
                    f"reaction:{reaction_id}"
                )
            precedents.append(
                {
                    "reaction_id": reaction_id,
                    "observation_id": observation_id,
                    "reaction_smiles": str(
                        context.get("reaction_smiles") or ""
                    ),
                    "reference_id": reference_id,
                    "reference_record": (
                        dict(reference) if reference is not None else None
                    ),
                    "experimental_detail": (
                        dict(experimental)
                        if experimental is not None
                        else None
                    ),
                }
            )
        recommendation["condition_precedents"] = precedents
    return payload


__all__ = ["attach_condition_precedents"]
