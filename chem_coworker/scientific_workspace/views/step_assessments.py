"""Project exact-structure saved assessments without rerunning scientific tools."""

from __future__ import annotations

from typing import Any

from ..adapters.step_selection import SOURCE_OPERATIONS, _selection
from ..answers.step_precedents import _structure_key
from ..core.store import InvestigationStore


def answer_step_assessments(
    store: InvestigationStore, answer: dict[str, Any], *, answer_ref: str | None = None,
) -> dict[str, dict[str, Any]]:
    """Expose all exact-step structural gates and cited recipe checks separately.

    Both complete sides, including stereo and chemical form, must match. Recipe
    checks describe their saved recipe, not necessarily the answer's prose recipe.
    Uncited structural checks cannot be hidden by citing only an inspection.
    An answer event bounds historical views so later checks do not rewrite them.
    Recipe checks remain citation-bound: matching graphs do not identify a recipe.
    """
    molecules = {item["id"]: item["smiles"] for item in answer.get("molecules", [])}
    result = {}
    keys = {}
    for step in answer.get("steps", []):
        keys[step["id"]] = _structure_key(
            ".".join(molecules[key] for key in step["reactant_ids"]),
            ".".join(molecules[key] for key in step["product_ids"]),
        )
        result[step["id"]] = {
            "schema_version": "scientific_step_evidence.v1",
            "status": "not_recorded", "structural_assessments": [], "recipe_assessments": [],
            "limitations": ["Saved assessments do not establish experimental feasibility."],
        }
    if not result:
        return result
    references = set(answer.get("evidence_refs", []))
    references.update(source["artifact_ref"] for source in answer.get("sources", []))
    references.update(ref for step in answer.get("steps", []) for ref in step.get("precedent_refs", []))
    references.update(ref for step in answer.get("steps", []) for ref in step.get("recipe_assessment_refs", []))
    try:
        events = store.events()
        if answer_ref is not None:
            boundary = next((event.sequence for event in events if event.artifact_ref == answer_ref), None)
            if boundary is None:
                raise ValueError("Saved answer event is unavailable")
            events = [event for event in events if event.sequence <= boundary]
        recorded_calls = {event.artifact_ref for event in events if event.kind == "call"}
        calls = {reference: store.read_artifact(reference)
                 for reference in sorted(references & recorded_calls)}
        # An attached inspection may cite its parent assessment indirectly.
        references.update(calls[ref]["result"]["source_ref"] for ref in tuple(references)
                          if ref in calls and calls[ref].get("operation") == "inspect_step_precedents")
        calls.update({reference: store.read_artifact(reference)
                      for reference in sorted((references & recorded_calls) - calls.keys())})
        # Recover assessment receipts independently of the agent's citation choices.
        # These are diagnostics, not newly inferred precedent or recipe support.
        for reference in sorted(recorded_calls - calls.keys()):
            payload = store.read_artifact(reference)
            if payload.get("operation") in SOURCE_OPERATIONS - {"disconnect_target"}:
                calls[reference] = payload
        for reference, payload in calls.items():
            if payload.get("execution_status") != "completed":
                continue
            operation = payload.get("operation")
            saved = payload.get("result", {})
            selections = []
            if operation in {"assess_recipe", "assess_proposed_recipe"}:
                sides = payload.get("arguments", {}).get("reaction_smiles", "").split(">")
                if len(sides) == 3:
                    selections.append((
                        {"precursor_smiles": sides[0], "target_smiles": sides[2]}, {},
                        saved.get("compatibility", saved),
                        saved.get("proposed_recipe", payload.get("arguments", {}).get("recipe", {})).get("recipe_id"),
                    ))
            elif operation in SOURCE_OPERATIONS and operation != "disconnect_target":
                selectors = [None] if operation in {"assess_route_step", "assess_retro_validity"} else [
                    item["external_step_id"] for item in saved.get("proposal", {}).get("steps", [])
                ]
                for selector in selectors:
                    selected, assessment = _selection(payload, selector, None)
                    recipe = (saved.get("proposed_recipe_assessment") if selector is None else
                              saved.get("proposed_recipe_assessments", {}).get(selector))
                    selections.append((selected, assessment, recipe, None))
            for selected, assessment, recipe, recipe_id in selections:
                key = _structure_key(selected["precursor_smiles"], selected["target_smiles"])
                if key is None:
                    continue
                for step_id, expected in keys.items():
                    if key != expected:
                        continue
                    view = result[step_id]
                    if assessment:
                        view["structural_assessments"].append({
                            "artifact_ref": reference,
                            "attribution": "cited" if reference in references else "saved_exact_step_check",
                            **{name: assessment.get(name) for name in (
                                "status", "actionable", "admission_eligible", "warnings",
                            )},
                            "gates": [{name: gate.get(name) for name in (
                                "gate_id", "status", "summary", "warnings",
                            )} for gate in assessment.get("gates", [])],
                        })
                    if recipe and reference in references:
                        view["recipe_assessments"].append({
                            "artifact_ref": reference, "recipe_id": recipe_id,
                            **({"process_coverage": saved.get("process_coverage")}
                               if operation == "assess_proposed_recipe" else {}),
                            **{name: recipe.get(name) for name in (
                                "status", "hard_conflicts", "checked_requirements",
                                "unresolved_requirements", "analysis_warnings", "coverage",
                                "evidence",
                            )},
                        })
                    if assessment or (recipe and reference in references):
                        view["status"] = "recorded"
    except (OSError, ValueError, KeyError, TypeError) as exc:
        # Fail closed rather than display a partial success after corrupt evidence.
        for view in result.values():
            view.update(status="evidence_unavailable", structural_assessments=[], recipe_assessments=[])
            view["limitations"].append(f"Saved assessment evidence could not be loaded: {exc}")
    return result
