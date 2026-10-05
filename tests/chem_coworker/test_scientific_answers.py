"""Versioned answer attribution and route-view integrity, not chemical feasibility tests."""

from __future__ import annotations

import base64
from copy import deepcopy
import json
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from app.web_api.scientific_presentation import present_conversation
from chem_coworker.scientific_workspace.answers.answer_contracts import (
    ANSWER_SCHEMA,
    ScientificAnswer,
    validate_answer_evidence,
)
from chem_coworker.scientific_workspace.core.store import InvestigationStore
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from chem_coworker.scientific_workspace.adapters.literature import capture_source
from chem_coworker.scientific_workspace.answers.evidence_review import answer_digest


def proposed() -> dict:
    return {"basis": "proposed", "source_ids": [], "limitations": ["Development illustration only"]}


@pytest.fixture
def example(tmp_path: Path):
    store = InvestigationStore.create(tmp_path / "investigation", objective="Presentation test", baseline={})
    excerpt = tmp_path / "user_excerpt.txt"
    excerpt.write_text(
        "Test fixture copied from user-provided text, not independently fetched or checked. "
        "Example 1: amino-group acylation; NaHCO3, THF, 16 h, 91.8%. "
        "Example 2: ring closure; p-TsOH, toluene, 16 h, 92.4%. "
        "The free-acid hydrolysis procedure was not included.", "utf-8",
    )
    source = store.attach_file(excerpt, description="User-quoted excerpt for a presentation regression")
    reported = {"basis": "reported", "source_ids": ["patent"], "limitations": ["User-quoted source; not independently verified"]}
    answer = {
        "schema_version": "scientific_answer.v2",
        "answer_markdown": "Illustrative structured view of the supplied excerpt; the route is incomplete.",
        "evidence_refs": [], "uncertainties": ["Hydrolysis procedure missing"], "needs_user_input": False,
        "sources": [{"id": "patent", "kind": "external_source", "title": "User-quoted patent examples",
                     "artifact_ref": source.artifact_ref, "url": "https://patents.google.com/patent/US12116352B2/en",
                     "locator": "User-supplied excerpt of Examples 1–2; not independently retrieved"}],
        "molecules": [
            {"id": "target", "name": "Target free acid", "smiles": "O=C(O)c1ccc2nc(-c3cc(Cl)cc(Cl)c3)oc2c1", "basis": "input", "source_ids": [], "limitations": []},
            {"id": "amino", "name": "Aminophenol ester", "smiles": "COC(=O)c1ccc(N)c(O)c1", **proposed()},
            {"id": "acid_chloride", "name": "Acid chloride", "smiles": "O=C(Cl)c1cc(Cl)cc(Cl)c1", **proposed()},
            {"id": "amide", "name": "Illustrative amide", "smiles": "COC(=O)c1ccc(NC(=O)c2cc(Cl)cc(Cl)c2)c(O)c1", **proposed()},
            {"id": "ester", "name": "Illustrative benzoxazole ester", "smiles": "COC(=O)c1ccc2nc(-c3cc(Cl)cc(Cl)c3)oc2c1", **proposed()},
        ],
        "target_molecule_ids": ["target"],
        "steps": [
            {"id": "s1", "title": "Acylation", "reactant_ids": ["amino", "acid_chloride"], "product_ids": ["amide"], "after_step_ids": [],
             "conditions": [{"text": "NaHCO3, dry THF, 16 h", **reported}], "yield_info": {"text": "91.8%", **reported}, **reported},
            {"id": "s2", "title": "Ring closure", "reactant_ids": ["amide"], "product_ids": ["ester"], "after_step_ids": ["s1"],
             "conditions": [{"text": "p-TsOH, toluene reflux, 16 h", **reported}], "yield_info": {"text": "92.4%", **reported}, **reported},
        ],
        "routes": [{"id": "r1", "title": "Illustration of the supplied two-step excerpt", "step_ids": ["s1", "s2"], "limitations": ["Does not yet reach the free acid"]}],
        "claims": [{"text": "Hydrolysis conditions are unknown from this excerpt", "basis": "unknown", "source_ids": ["patent"], "limitations": []}],
    }
    return store, answer


def test_tafamidis_view_renders_steps_and_retains_incomplete_route(example) -> None:
    store, payload = example
    answer = ScientificAnswer.model_validate(payload)
    assert validate_answer_evidence(answer, store) == [payload["sources"][0]["artifact_ref"]]
    conversation = {"id": "a" * 32, "turns": [{"question": "How to make this target?", "answer": payload}]}
    original = deepcopy(conversation)
    view = present_conversation(conversation)["turns"][0]["structured_presentation"]
    assert conversation == original
    assert view["routes"][0]["unreached_target_ids"] == ["target"]
    assert view["routes"][0]["drawing_status"] == "linear_scheme"
    route_svg = ET.fromstring(base64.b64decode(view["routes"][0]["image_url"].split(",", 1)[1]))
    assert route_svg.get("data-layout") == "three-column-row-major"
    assert all(step["drawing_status"] == "drawn" for step in view["steps"])
    assert all(molecule["drawing_status"] == "drawn" for molecule in view["molecules"])
    for item in view["routes"] + view["steps"]:
        root = ET.fromstring(base64.b64decode(item["image_url"].split(",", 1)[1]))
        assert root.tag.endswith("svg")
    assert view["steps"][0]["yield_info"]["basis"] == "reported"
    assert view["molecules"][3]["basis"] == "proposed"
    scheme = ET.fromstring(base64.b64decode(view["steps"][0]["image_url"].split(",", 1)[1]))
    text = " ".join(scheme.itertext())
    assert "NaHCO3, dry THF, 16 h" in text and "Yield (reported): 91.8%" in text
    assert "Aminophenol ester" in text and "Illustrative amide" in text


@pytest.mark.parametrize("case", ["branch", "invalid", "ambiguous_product"])
def test_route_scheme_fallback_preserves_steps_and_explains_reason(example, case) -> None:
    _, payload = example
    if case == "branch":
        branch = deepcopy(payload["steps"][1])
        branch["id"] = "s3"
        payload["steps"].append(branch)
        payload["routes"][0]["step_ids"].append("s3")
    elif case == "invalid":
        payload["molecules"][1]["smiles"] = "C1"
    else:
        payload["steps"][1]["product_ids"].append("target")
    original = deepcopy(payload)
    view = present_conversation({"id": "a" * 32, "turns": [
        {"question": "Q", "answer": payload},
    ]})["turns"][0]["structured_presentation"]
    route = view["routes"][0]
    assert route["drawing_status"] == "dependency_overview"
    assert route["drawing_warning"]
    assert "scheme_width" not in route
    assert len(view["steps"]) == len(payload["steps"])
    assert payload == original


def test_saved_answer_uses_source_title_and_real_url_without_mutation(example) -> None:
    _, payload = example
    reference = payload["sources"][0]["artifact_ref"]
    payload["answer_markdown"] = f"See {reference} for the procedure."
    original = deepcopy(payload)
    view = present_conversation({"id": "a" * 32, "turns": [{"question": "Q", "answer": payload}]})
    html = view["turns"][0]["answer_presentation"]["html"]
    assert "User-quoted patent examples</a>" in html
    assert 'href="https://patents.google.com/patent/US12116352B2/en"' in html
    assert payload == original


def test_route_labels_use_explicit_reagents_and_keep_full_conditions(example) -> None:
    _, payload = example
    step = payload["steps"][0]
    step["reagents"] = [{**step["conditions"][0], "text": "NaHCO3"}]
    original = deepcopy(payload)
    view = present_conversation({"id": "a" * 32, "turns": [
        {"question": "Q", "answer": payload},
    ]})["turns"][0]["structured_presentation"]
    root = ET.fromstring(base64.b64decode(view["routes"][0]["image_url"].split(",", 1)[1]))
    ns = "{http://www.w3.org/2000/svg}"
    visible = " ".join(node.text or "" for node in root.iter(ns + "text"))
    assert "NaHCO3" in visible
    assert "THF" not in visible and "16 h" not in visible and "91.8%" not in visible
    assert "p-TsOH" not in visible  # No reagent names guessed from old procedure text.
    assert view["steps"][0]["conditions"] == step["conditions"]
    assert view["steps"][0]["yield_info"] == step["yield_info"]
    assert payload == original


@pytest.mark.parametrize("source_ids,match", [([], "require source_ids"), (["absent"], "Unknown source_id")])
def test_reagent_claims_require_valid_attribution(example, source_ids, match) -> None:
    _, payload = example
    payload["steps"][0]["reagents"] = [{
        "text": "NaHCO3", "basis": "reported", "source_ids": source_ids, "limitations": [],
    }]
    with pytest.raises(ValueError, match=match):
        ScientificAnswer.model_validate(payload)


@pytest.mark.parametrize("mutation,match", [
    (lambda a: a["steps"][0].update(source_ids=[]), "require source_ids"),
    (lambda a: a["steps"][0].update(source_ids=["absent"]), "Unknown source_id"),
    (lambda a: a["steps"][0].update(product_ids=["absent"]), "Unknown molecule"),
    (lambda a: a["steps"][0].update(after_step_ids=["s2"]), "acyclic"),
    (lambda a: a["steps"][1].update(reactant_ids=["target"]), "intermediate"),
    (lambda a: a["routes"][0].update(step_ids=["s2"]), "preceding dependencies"),
    (lambda a: a["molecules"][1].update(id="target"), "unique"),
    (lambda a: a["sources"][0].update(url="javascript:alert(1)"), "HTTP"),
    (lambda a: a["steps"][0].update(basis="computed"), "local computation"),
])
def test_bad_attribution_or_route_connections_are_rejected(example, mutation, match) -> None:
    _, payload = example
    mutation(payload)
    with pytest.raises(ValueError, match=match):
        ScientificAnswer.model_validate(payload)


def test_external_reference_cannot_be_just_a_url_or_an_agent_note(example) -> None:
    store, payload = example
    note = store.note("hypothesis", "This is not a captured source")
    payload["sources"][0]["artifact_ref"] = note.artifact_ref
    with pytest.raises(ValueError, match="scientific evidence"):
        validate_answer_evidence(ScientificAnswer.model_validate(payload), store)


def test_computation_requires_completed_call_but_does_not_override_domain_validity(example) -> None:
    store, payload = example
    call = store.append("call", {"operation": "analyze_reaction", "execution_status": "error", "error": "unsupported"})
    payload["sources"] = [{"id": "local", "kind": "local_artifact", "title": "Recorded analysis", "artifact_ref": call.artifact_ref, "url": None, "locator": "result.valid"}]
    payload.update(molecules=[], target_molecule_ids=[], steps=[], routes=[], claims=[{
        "text": "No supported structural interpretation", "basis": "computed", "source_ids": ["local"], "limitations": ["Not evidence of chemical impossibility"],
    }])
    with pytest.raises(ValueError, match="completed recorded"):
        validate_answer_evidence(ScientificAnswer.model_validate(payload), store)
    attachment = store.append("derived_file", {
        "operation": "analyze_reaction", "execution_status": "completed", "result": {"valid": True},
    })
    payload["sources"][0]["artifact_ref"] = attachment.artifact_ref
    with pytest.raises(ValueError, match="completed recorded"):
        validate_answer_evidence(ScientificAnswer.model_validate(payload), store)
    completed = store.append("call", {"operation": "analyze_reaction", "execution_status": "completed", "result": {"valid": False}})
    payload["sources"][0]["artifact_ref"] = completed.artifact_ref
    assert validate_answer_evidence(ScientificAnswer.model_validate(payload), store) == [completed.artifact_ref]


def test_invalid_smiles_remains_inspectable_and_old_answers_are_not_reconstructed(example) -> None:
    _, payload = example
    payload["molecules"][1]["smiles"] = "INVALID"
    view = present_conversation({"id": "a" * 32, "turns": [{"question": "Q", "answer": payload}]})
    structured = view["turns"][0]["structured_presentation"]
    assert structured["molecules"][1]["drawing_status"] == "invalid_or_unsupported_notation"
    assert structured["steps"][0]["drawing_status"] == "invalid_or_unsupported_notation"
    legacy = present_conversation({"id": "a" * 32, "turns": [{"question": "Q", "answer": {"answer_markdown": "Old saved prose"}}]})
    assert "structured_presentation" not in legacy["turns"][0]


def test_runtime_schema_requires_closed_explicit_objects_and_nulls() -> None:
    for schema in [ANSWER_SCHEMA, *ANSWER_SCHEMA["$defs"].values()]:
        if schema.get("type") == "object":
            assert schema["additionalProperties"] is False
            assert set(schema["required"]) == set(schema["properties"])
    assert json.loads(json.dumps(ANSWER_SCHEMA))["properties"]["schema_version"]["const"] == "scientific_answer.v2"


def literature_example(example):
    """Attach a distinct source reaction, rather than copy the target proposal."""
    store, payload = example
    source = payload["sources"][0]
    event = capture_source(store, "Example 9: ethanol was oxidized to acetaldehyde. CCO>>CC=O.", url=source["url"])
    source["artifact_ref"] = event.artifact_ref
    drawing = {
        "title": "Source oxidation", "source_id": "patent", "locator": "Example 9",
        "structure_origin": "reconstructed_from_description",
        "structure_evidence": "ethanol was oxidized to acetaldehyde.",
        "reactants": [{"name": "Ethanol", "smiles": "CCO"}],
        "products": [{"name": "Acetaldehyde", "smiles": "CC=O"}],
        "limitations": ["Source structure interpretation remains unverified."],
    }
    payload["steps"][0]["literature_reactions"] = [drawing]
    return store, payload, drawing


@pytest.mark.parametrize("origin", ["source_explicit", "reconstructed_from_description"])
def test_source_reactions_render_automatically_without_copying_target(example, origin):
    store, payload, drawing = literature_example(example)
    drawing["structure_origin"] = origin
    validate_answer_evidence(ScientificAnswer.model_validate(payload), store)
    original = deepcopy(payload)
    view = present_conversation({"id": "a" * 32, "turns": [{"question": "Q", "answer": payload}]})["turns"][0]["structured_presentation"]
    source = view["steps"][0]["literature_reactions"][0]
    assert source["reaction_smiles"] == "CCO>>CC=O"
    assert source["drawing_status"] == "drawn"
    assert source["reaction_smiles"] != view["steps"][0]["reaction_smiles"]
    assert source["limitations"] == drawing["limitations"]
    ET.fromstring(base64.b64decode(source["image_url"].split(",", 1)[1]))
    assert payload == original


@pytest.mark.parametrize("change,match", [
    ({"source_id": "missing"}, "captured external"),
    ({"structure_evidence": "Invented passage"}, "quote its captured source"),
    ({"structure_evidence": "  "}, "must not be blank"),
    ({"structure_origin": "source_explicit", "products": [{"name": "Different product", "smiles": "CCN"}]}, "SMILES must appear"),
    ({"conditions": [{"text": "Target adaptation", "basis": "proposed", "source_ids": ["patent"], "limitations": []}]}, "must be reported"),
])
def test_literature_drawing_rejects_missing_and_conflicting_provenance(example, change, match):
    store, payload, drawing = literature_example(example)
    drawing.update(change)
    with pytest.raises(ValueError, match=match):
        validate_answer_evidence(ScientificAnswer.model_validate(payload), store)


def test_invalid_source_structure_is_retained_without_a_drawing(example):
    _, payload, drawing = literature_example(example)
    drawing["products"][0]["smiles"] = "INVALID"
    view = present_conversation({"id": "a" * 32, "turns": [{"question": "Q", "answer": payload}]})["turns"][0]["structured_presentation"]
    source = view["steps"][0]["literature_reactions"][0]
    assert source["drawing_status"] == "invalid_or_unsupported_notation"
    assert "image_url" not in source
    assert source["products"][0]["smiles"] == "INVALID"


def test_optional_source_drawings_preserve_old_review_hashes_and_bind_new_drawings(example):
    _, payload = example
    old = ScientificAnswer.model_validate(payload).model_dump()
    old["evidence_refs"] = sorted(set(old["evidence_refs"]))
    for step in old["steps"]:
        step.pop("literature_reactions")
        # Historical drafts predate both optional evidence extensions.
        step.pop("recipe_assessment_refs")
    assert answer_digest(payload) == hashlib.sha256(canonical_bytes(old)).hexdigest()
    _, payload, drawing = literature_example(example)
    initial = answer_digest(payload)
    drawing["products"][0]["smiles"] = "CCN"
    assert answer_digest(payload) != initial
