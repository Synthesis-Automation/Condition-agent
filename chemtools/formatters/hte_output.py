"""
HTE recommendation output formatting.

Builds the concise, execution-focused JSON structure defined by
chemtools/formatters/output_core_format.json.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from chemtools import reagent


_ROLE_TO_REAGENT_TYPE = {
    "metal_catalyst": "metal_catalyst",
    "catalyst": "metal_catalyst",
    "ligand": "ligand",
    "base": "base",
    "solvent": "solvent",
    "additive": "additive",
    "coupling_reagent": "condensation_agent",
}

_DEFAULT_EQUIVALENTS_BY_ROLE = {
    "metal_catalyst": 0.05,
    "catalyst": 0.05,
    "ligand": 0.05,
    "base": 2.0,
    "acid": 2.0,
    "additive": 1.5,
    "coupling_reagent": 1.5,
}
_DEFAULT_OTHER_EQUIVALENTS = 1.5
_DEFAULT_SOLVENT_CONCENTRATION_M = 0.2


def _clean_text(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _reaction_smiles_from_result(result: Any) -> str:
    reactants: List[str] = []
    reactant_a = _clean_text(getattr(result, "reactant_a_smiles", None))
    reactant_b = _clean_text(getattr(result, "reactant_b_smiles", None))
    product = _clean_text(getattr(result, "product_smiles", None))
    if reactant_a:
        reactants.append(reactant_a)
    if reactant_b:
        reactants.append(reactant_b)
    reactants_text = ".".join(reactants)
    if product:
        return f"{reactants_text}>>{product}" if reactants_text else f">>{product}"
    return reactants_text


def _lookup_reagent_info(name: Optional[str], role: str) -> Dict[str, Any]:
    cleaned = _clean_text(name)
    if not cleaned:
        return {"name": None, "cas": None, "smiles": None}
    reagent_type = _ROLE_TO_REAGENT_TYPE.get(role, role)
    try:
        info = reagent.enrich_reagent_info(cleaned, reagent_type)
    except Exception:
        info = {}
    if info.get("found"):
        return {
            "name": info.get("name") or cleaned,
            "cas": info.get("cas"),
            "smiles": info.get("smiles"),
        }
    return {"name": cleaned, "cas": None, "smiles": None}


def _chemical_entry(
    *,
    name: Optional[str],
    role: str,
    cas: Optional[str] = None,
    smiles: Optional[str] = None,
    equivalents: Optional[float] = None,
    extras: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "name": name,
        "cas": cas,
        "smiles": smiles,
        "equivalents": equivalents,
        "role": role,
    }
    if extras:
        entry.update(extras)
    if entry.get("equivalents") is None:
        is_limiting = bool(entry.get("is_limiting_reactant"))
        if role == "reactant":
            entry["equivalents"] = 1.0 if is_limiting else 1.5
        elif role != "solvent":
            entry["equivalents"] = _DEFAULT_EQUIVALENTS_BY_ROLE.get(
                role, _DEFAULT_OTHER_EQUIVALENTS
            )
    if role == "solvent" and entry.get("concentration_M") is None:
        entry["concentration_M"] = _DEFAULT_SOLVENT_CONCENTRATION_M
    return entry


def _append_add_step(steps: List[Dict[str, Any]], chemical: Dict[str, Any]) -> None:
    if not chemical.get("name") and not chemical.get("smiles"):
        return
    steps.append({"action": "add", "chemical": chemical})


def format_hte_output(
    result: Any,
    *,
    recommendations: Optional[List[Any]] = None,
    reaction_smiles: Optional[str] = None,
    reaction_type_filter: Optional[str] = None,
    catalyst_filter: Optional[str] = None,
    explanation: Optional[str] = None,
    condition_steps: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """
    Convert HTERecommendationResult into the concise execution JSON format.

    Args:
        result: HTERecommendationResult instance.
        reaction_smiles: Optional override for the reaction SMILES.
        reaction_type_filter: Optional user-specified reaction filter.
        catalyst_filter: Optional user-specified catalyst filter.
        explanation: Optional explanation text for the recommendation.
        condition_steps: Optional list of condition step dicts to append.

    Returns:
        JSON-ready dictionary matching output_core_format.json structure.
    """
    reaction_smiles = reaction_smiles or _reaction_smiles_from_result(result)
    input_block: Dict[str, Any] = {"reaction_smiles": reaction_smiles}
    if reaction_type_filter:
        input_block["reaction_type_filter"] = reaction_type_filter
    if catalyst_filter:
        input_block["catalyst_filter"] = catalyst_filter

    detection_block: Dict[str, Any] = {}
    predicted = _clean_text(getattr(result, "predicted_reaction_type", None))
    if predicted:
        detection_block["reaction_type"] = predicted
    confidence = getattr(result, "reaction_type_confidence", None)
    if isinstance(confidence, (int, float)):
        detection_block["confidence"] = round(float(confidence), 4)

    output: Dict[str, Any] = {
        "meta": {},
        "input": input_block,
        "detection": detection_block,
        "explanation": {"text": explanation or ""},
        "recommended_conditions": [],
    }

    reactant_a_smiles = _clean_text(getattr(result, "reactant_a_smiles", None))
    reactant_b_smiles = _clean_text(getattr(result, "reactant_b_smiles", None))

    recs = list(recommendations or (getattr(result, "recommendations", []) or []))
    for index, rec in enumerate(recs, start=1):
        steps: List[Dict[str, Any]] = []

        if reactant_a_smiles:
            _append_add_step(
                steps,
                _chemical_entry(
                    name=None,
                    cas=None,
                    smiles=reactant_a_smiles,
                    equivalents=None,
                    role="reactant",
                    extras={"is_limiting_reactant": True},
                ),
            )

        if reactant_b_smiles:
            _append_add_step(
                steps,
                _chemical_entry(
                    name=None,
                    cas=None,
                    smiles=reactant_b_smiles,
                    equivalents=None,
                    role="reactant",
                ),
            )

        solvent = _clean_text(getattr(rec, "solvent", None))
        if solvent:
            info = _lookup_reagent_info(solvent, "solvent")
            _append_add_step(
                steps,
                _chemical_entry(
                    name=info["name"],
                    cas=info["cas"],
                    smiles=info["smiles"],
                    equivalents=None,
                    role="solvent",
                ),
            )

        secondary = _clean_text(getattr(rec, "secondary_solvent", None))
        if secondary:
            info = _lookup_reagent_info(secondary, "solvent")
            _append_add_step(
                steps,
                _chemical_entry(
                    name=info["name"],
                    cas=info["cas"],
                    smiles=info["smiles"],
                    equivalents=None,
                    role="solvent",
                    extras={"is_secondary_solvent": True},
                ),
            )

        base = _clean_text(getattr(rec, "base", None))
        if base:
            info = _lookup_reagent_info(base, "base")
            _append_add_step(
                steps,
                _chemical_entry(
                    name=info["name"],
                    cas=info["cas"],
                    smiles=info["smiles"],
                    equivalents=None,
                    role="base",
                ),
            )

        catalyst = _clean_text(getattr(rec, "catalyst", None))
        if catalyst:
            info = _lookup_reagent_info(catalyst, "metal_catalyst")
            _append_add_step(
                steps,
                _chemical_entry(
                    name=info["name"],
                    cas=info["cas"],
                    smiles=info["smiles"],
                    equivalents=None,
                    role="metal_catalyst",
                ),
            )

        ligand = _clean_text(getattr(rec, "ligand", None))
        if ligand:
            info = _lookup_reagent_info(ligand, "ligand")
            _append_add_step(
                steps,
                _chemical_entry(
                    name=info["name"],
                    cas=info["cas"],
                    smiles=info["smiles"],
                    equivalents=None,
                    role="ligand",
                ),
            )

        additive = _clean_text(getattr(rec, "additive", None))
        if additive:
            info = _lookup_reagent_info(additive, "additive")
            _append_add_step(
                steps,
                _chemical_entry(
                    name=info["name"],
                    cas=info["cas"],
                    smiles=info["smiles"],
                    equivalents=None,
                    role="additive",
                ),
            )

        coupling_reagent = _clean_text(getattr(rec, "coupling_reagent", None))
        if coupling_reagent:
            info = _lookup_reagent_info(coupling_reagent, "coupling_reagent")
            _append_add_step(
                steps,
                _chemical_entry(
                    name=info["name"],
                    cas=info["cas"],
                    smiles=info["smiles"],
                    equivalents=None,
                    role="coupling_reagent",
                ),
            )

        if condition_steps:
            for condition in condition_steps:
                if not isinstance(condition, dict):
                    continue
                step = {"action": "condition"}
                step.update(condition)
                steps.append(step)

        output["recommended_conditions"].append({"rank": index, "steps": steps})

    return output
