"""Explicit public helper documentation for agents; no implementation-file discovery."""

from __future__ import annotations

from typing import Any

HELPER_EXAMPLES = {
    "run": "event = w.run(operation_name, parameters)",
    "run_summary": "print(w.run_summary(operation_name, parameters))",
    "call_summary": "print(w.call_summary(saved_ref))",
    "fetch_source": "source = w.fetch_source(paper_url); print(w.call_summary(source))",
    "capture_source_file": "source = w.capture_source_file('browser.json', url=paper_url, text_path=['text'])",
    "capture_source": "source = w.capture_source(text, url=paper_url, locator='Example 1')",
    "capture_source_image": "image = w.capture_source_image('scheme.png', source_ref=source_ref, locator='Scheme 1, page 2')",
    "record_source_excerpt": "excerpt = w.record_source_excerpt(source_ref, start=120, end=640, locator='Example 1')",
    "inspect_source": "print(w.inspect_source(source_ref, query='Compound 8', limit=2000))",
    "batch_summary": "print(w.batch_summary([event.artifact_ref for event in events]))",
    "prepared_literature_reaction": "reaction = w.prepared_literature_reaction(preparation_ref)",
    "attach_literature_reaction": "draft = w.attach_literature_reaction(draft, 's1', preparation_ref)",
    "attach_recipe_check": "draft = w.attach_recipe_check(draft, 's1', recipe_check_ref)",
    "answer_template": "draft = w.answer_template('An evidence-linked proposal with explicit gaps.')",
    "route_answer": "draft = w.route_answer({'target_smiles': 'CC=O', 'routes': [{'steps': [{'reaction_smiles': 'CCO>>CC=O'}]}]})",
    "inspect_artifact": "print(w.inspect_artifact(ref, path=('result',), offset=0, limit=3))",
    "answer_preflight": "print(w.answer_preflight(draft))",
    "finalize_answer": "receipt = w.finalize_answer(draft_path, draft, findings=findings)",
    "run_python": "event = w.run_python('audit.py', parameters, evidence_refs=(source_ref,))",
}


def nested_input_help(name: str) -> dict[str, Any]:
    """Expose actual nested contracts only when their operation is requested."""
    if name == "route_answer":
        from ..answers.route_authoring import RouteAnswerInput

        return {"input_schema": RouteAnswerInput.model_json_schema(), "input_notes": [
            "Supply reactants>>products and inspected support as ref/locator; conditions, titles and notes are optional.",
            "Steps and conditions default to proposed. Reported labels require support; missing conditions stay missing.",
            "after_steps optionally names earlier steps by one-based route-local position; otherwise unique exact intermediate links are assembled.",
            "No science or source drawings run. finalize_answer still checks citations and required exact-step inspections.",
        ]}
    if name == "assess_proposed_recipe":
        from pydantic import TypeAdapter
        from condition_registry.models import ConditionComponentInput, ConditionProcessStage

        return {
            "nested_inputs": {
                "components[]": TypeAdapter(ConditionComponentInput).json_schema(),
                "operating_conditions.stages[]": TypeAdapter(ConditionProcessStage).json_schema(),
            },
            "example_arguments": {
                "reaction_smiles": "CCBr.N>>CCN",
                "components": [{"raw_identifier": "ethanol", "source_field": "proposal",
                                "source_role_hint": "solvent", "provenance": {"description": "Proposed medium"}}],
                "operating_conditions": {},
            },
            "input_notes": ["provenance is a JSON object, not text.",
                            "Represent all evidence-backed atom contributors and required multiplicity in reaction_smiles. "
                            "Recipe amounts do not add atoms to that graph; never invent a donor or mapping.",
                            "Try one input before batching the same new nested format."],
        }
    if name == "prepare_literature_reaction":
        from ..adapters.literature_reactions import ParticipantInput, SourceConflict
        from ..answers.answer_contracts import AnswerClaim

        return {"nested_inputs": {"reactants[] / products[]": ParticipantInput.model_json_schema(),
                                  "conditions[] / yield_info": AnswerClaim.model_json_schema(),
                                  "source_conflicts[]": SourceConflict.model_json_schema()},
                "input_notes": ["conditions and yield_info require reported basis, source_ids and literal captured text; "
                                "omitted limitations default to an empty list.",
                                "Adjacent mass/amount pairs are checked conditionally against supplied graphs. "
                                "Empty quantity_checks means no eligible pair was found, not consistency."]}
    if name in {"answer_template", "answer_preflight", "finalize_answer"}:
        from ..answers.answer_contracts import AnswerClaim

        return {"nested_inputs": {"claims[]": AnswerClaim.model_json_schema()},
                "input_notes": ["Claims have text, basis, source_ids and limitations; no id field.",
                                "Correct a failed answer_preflight before calling finalize_answer; "
                                "a valid draft can still have unresolved scientific warnings."]}
    return {}
