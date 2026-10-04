"""Exact participant evidence, indexed reuse, conflicts and immutable answer preparation."""

from copy import deepcopy
import base64
import json
from pathlib import Path

import pytest
from PIL import Image

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.answers.answer_contracts import ScientificAnswer, validate_answer_evidence
from chem_coworker.scientific_workspace.core.baseline import code_manifest, environment_versions
from chem_coworker.scientific_workspace.core.store import canonical_bytes
from chem_coworker.scientific_workspace.__main__ import main
from chem_coworker.scientific_workspace.runtime.activity import ActivityHistory, ScientificActivityCursor
from app.web_api.scientific_presentation import present_conversation

ROOT = Path(__file__).resolve().parents[2]
PUBLICATION = "REF1:" + "a" * 64
TEXT = "Compound 8: bromoethane C2H5Br; ammonia N. Compound 9: ethylamine C2H7N.\nDiscussion: DCM.\nExperimental: DCE."


def test_source_image_is_preserved_bound_and_presented_without_verifying_assignment(workspace):
    from chem_coworker.scientific_workspace.views.literature_images import answer_literature_images

    source, args = setup(workspace)
    path = workspace.store.root / "scheme.png"
    Image.new("RGB", (10, 10), "white").save(path)
    original = path.read_bytes()
    image = workspace.capture_source_image("scheme.png", source_ref=source.artifact_ref, locator="Scheme 1")
    assert "image_base64" not in str(workspace.call_summary(image))
    args["scheme_refs"] = [image.artifact_ref]
    preparation, block = prepare(workspace, args)
    draft = answer(source, block).model_dump()
    evidence = answer_literature_images(workspace.store, draft)
    captured = evidence[preparation.artifact_ref]["images"][0]
    assert base64.b64decode(captured["image_url"].split(",", 1)[1]) == original
    assert captured["assignment_verification"] == "not_performed"
    before = deepcopy(draft)
    shown = present_conversation({"id": "a" * 32, "turns": [{"id": "t", "question": "prepare",
        "answer": draft, "literature_image_evidence": evidence}]})["turns"][0]["structured_presentation"]
    assert shown["steps"][0]["literature_reactions"][0]["captured_images"][0] == captured
    assert draft == before
    other = workspace.capture_source("Other source", url="https://example.org/other")
    image2 = workspace.capture_source_image("scheme.png", source_ref=other.artifact_ref, locator="Other scheme")
    args["scheme_refs"] = [image2.artifact_ref]
    failed = workspace.run("prepare_literature_reaction", args)
    assert workspace.store.read_artifact(failed.artifact_ref)["execution_status"] == "error"


def test_attachment_helper_uses_exact_preparation_source_and_preserves_draft(workspace):
    source, args = setup(workspace)
    event, block = prepare(workspace, args)
    draft = answer(source, block).model_dump()
    draft["steps"][0]["literature_reactions"] = []
    draft["sources"] = []
    before = deepcopy(draft)
    attached = workspace.attach_literature_reaction(draft, "s", event.artifact_ref)
    assert attached["sources"][0]["artifact_ref"] == source.artifact_ref
    validate_answer_evidence(ScientificAnswer.model_validate(attached), workspace.store)
    assert draft == before
    conflicting = deepcopy(attached)
    conflicting["sources"][0]["artifact_ref"] = args["reactants"][0]["evidence_ref"]
    with pytest.raises(ValueError, match="conflicts"):
        workspace.attach_literature_reaction(conflicting, "s", event.artifact_ref)


@pytest.fixture
def workspace(tmp_path: Path) -> ScientificWorkspace:
    store = InvestigationStore.create(tmp_path / "study", objective="Prepare source reaction", baseline={
        "repository": str(ROOT), "code_files": code_manifest(ROOT),
        "environment": environment_versions(), "artifacts": {},
    })
    return ScientificWorkspace(store.root)


def setup(workspace: ScientificWorkspace):
    source = workspace.capture_source(TEXT, url="https://example.org/paper", reference_id=PUBLICATION)
    refs = [workspace.record_source_excerpt(source.artifact_ref, excerpt=quote).artifact_ref for quote in
            ["Compound 8: bromoethane C2H5Br; ammonia N.", "Compound 9: ethylamine C2H7N."]]
    arguments = {
        "source_ref": source.artifact_ref, "source_id": "paper", "title": "Source reaction",
        "locator": "Example 1", "structure_evidence": "Compound 9: ethylamine C2H7N.",
        "reactants": [
            {"name": "Bromoethane", "compound_id": "8", "evidence_ref": refs[0],
             "material_form": "Neutral molecule", "smiles": "CCBr", "reported_formula": "C2H5Br"},
            {"name": "Ammonia", "compound_id": "ammonia", "evidence_ref": refs[0],
             "material_form": "Source concentration unavailable", "smiles": "N"},
        ],
        "products": [{"name": "Ethylamine", "compound_id": "9", "evidence_ref": refs[1],
                      "material_form": "Free base; isolation form unverified", "smiles": "CCN",
                      "reported_formula": "C2H7N"}],
    }
    return source, arguments


def prepare(workspace: ScientificWorkspace, arguments: dict):
    event = workspace.run("prepare_literature_reaction", arguments)
    record = workspace.store.read_artifact(event.artifact_ref)
    assert record["execution_status"] == "completed", record.get("error")
    return event, workspace.prepared_literature_reaction(event.artifact_ref)


def answer(source, block):
    return ScientificAnswer.model_validate({
        "schema_version": "scientific_answer.v2", "answer_markdown": "Proposed route; source assignment unverified.",
        "evidence_refs": [], "uncertainties": ["Source assignment unverified"], "needs_user_input": False,
        "sources": [{"id": "paper", "kind": "external_source", "title": "Captured paper",
                     "artifact_ref": source.artifact_ref, "url": "https://example.org/paper", "locator": "Example 1"}],
        "molecules": [{"id": "r", "name": "Bromoethane", "smiles": "CCBr", "basis": "proposed",
                       "source_ids": [], "limitations": []},
                      {"id": "p", "name": "Ethylamine", "smiles": "CCN", "basis": "proposed",
                       "source_ids": [], "limitations": []}],
        "target_molecule_ids": ["p"], "steps": [{"id": "s", "title": "Proposed step",
            "reactant_ids": ["r"], "product_ids": ["p"], "after_step_ids": [], "conditions": [],
            "yield_info": None, "basis": "proposed", "source_ids": [], "limitations": [],
            "literature_reactions": [block]}], "routes": [], "claims": [],
    })


def test_preparation_preserves_exact_evidence_and_replays(workspace):
    source, arguments = setup(workspace)
    event, block = prepare(workspace, arguments)
    assert block["schema_version"] == "literature_reaction.v2"
    assert block["source_provenance"]["acquisition"] == "agent_supplied_excerpt"
    assert block["products"][0]["formula"] == "C2H7N"
    assert block["products"][0]["graph_status"] == "graph_checked_assignment_unverified"
    assert arguments["products"][0]["evidence_ref"] in event.evidence_refs
    validate_answer_evidence(answer(source, block), workspace.store)
    replay = workspace.replay(event.artifact_ref)
    assert workspace.store.read_artifact(replay.artifact_ref)["matches"]


def test_source_quantity_conflict_survives_preparation_attachment_and_preflight(workspace):
    text = "Methyl iodide (162 g, 2.60 mol). Product P: ethane."
    source = workspace.capture_source(text, url="https://example.org/quantity")
    excerpt = workspace.record_source_excerpt(source.artifact_ref, excerpt=text)
    event, block = prepare(workspace, {
        "source_ref": source.artifact_ref, "source_id": "paper", "title": "Quantity discrepancy",
        "locator": "Example", "structure_evidence": "Product P: ethane.",
        "reactants": [{"name": "Methyl iodide", "compound_id": "Methyl iodide", "smiles": "CI",
                       "evidence_ref": excerpt.artifact_ref, "material_form": "Named reagent; assignment unverified"}],
        "products": [{"name": "Ethane", "compound_id": "P", "smiles": "CC",
                      "evidence_ref": excerpt.artifact_ref, "material_form": "Named product"}],
        "conditions": [{"text": "Methyl iodide (162 g, 2.60 mol).", "basis": "reported", "source_ids": ["paper"]}],
    })
    saved = workspace.store.read_artifact(event.artifact_ref)["result"]
    assert saved["schema_version"] == "literature_reaction_preparation.v2"
    assert saved["quantity_checks"][0]["status"] == "conflicting"
    assert saved["quantity_checks"][0]["evidence_ref"] == excerpt.artifact_ref
    assert any("mass/amount conflict" in note for note in block["limitations"])
    assert block["conditions"][0]["limitations"] == []
    brief = workspace.call_summary(event)["result_summary"]
    assert brief["quantity_conflict_count"] == 1
    draft = answer(source, block).model_dump()
    draft["sources"][0]["url"] = "https://example.org/quantity"
    preflight = workspace.answer_preflight(draft)
    assert preflight["valid"]
    assert any(item["gap"] == "source_quantity_conflict" for item in preflight["warnings"])
    assert workspace.store.read_artifact(workspace.replay(event.artifact_ref).artifact_ref)["matches"]


def test_help_exposes_existing_fetch_helper_and_actual_nested_contracts(workspace):
    entries = {item["name"]: item for item in workspace.help([
        "fetch_source", "assess_proposed_recipe", "prepare_literature_reaction", "finalize_answer",
    ])}
    assert "retry_network" in entries["fetch_source"]["signature"]
    schema = entries["assess_proposed_recipe"]["nested_inputs"]["components[]"]
    assert schema["properties"]["provenance"]["type"] == "object"
    example = entries["assess_proposed_recipe"]["example_arguments"]
    event = workspace.run("assess_proposed_recipe", example)
    assert workspace.store.read_artifact(event.artifact_ref)["execution_status"] == "completed"
    assert "compound_id" in entries["prepare_literature_reaction"]["nested_inputs"]["reactants[] / products[]"]["required"]
    claims = entries["finalize_answer"]["nested_inputs"]["claims[]"]
    assert "id" not in claims["properties"] and claims["additionalProperties"] is False


@pytest.mark.parametrize("case", ["compound", "formula", "other_source", "fake_explicit", "missing_evidence"])
def test_unsupported_attribution_is_rejected(workspace, case):
    _, arguments = setup(workspace)
    if case == "compound":
        arguments["products"][0]["compound_id"] = "UNKNOWN_COMPOUND"
    elif case == "formula":
        arguments["products"][0]["reported_formula"] = "C99H99N"
    elif case == "other_source":
        other = workspace.capture_source("Compound 9", url="https://example.org/other")
        arguments["products"][0]["evidence_ref"] = workspace.record_source_excerpt(
            other.artifact_ref, excerpt="Compound 9").artifact_ref
    elif case == "fake_explicit":
        arguments["structure_origin"] = "source_explicit"
    else:
        arguments["products"][0].pop("evidence_ref")
    record = workspace.store.read_artifact(workspace.run("prepare_literature_reaction", arguments).artifact_ref)
    assert record["execution_status"] == "error"


@pytest.mark.parametrize("smiles,status", [("CCC", "conflicting"), ("C1", "invalid")])
def test_formula_conflict_and_invalid_graph_are_not_repaired(workspace, smiles, status):
    _, arguments = setup(workspace)
    arguments["products"][0]["smiles"] = smiles
    _, block = prepare(workspace, arguments)
    assert block["products"][0]["smiles"] == smiles
    assert block["products"][0]["graph_status"] == status
    assert any("formula conflict" in note for note in block["limitations"])


def test_source_discrepancy_retains_both_passages(workspace):
    source, arguments = setup(workspace)
    refs = [workspace.record_source_excerpt(source.artifact_ref, excerpt=quote).artifact_ref
            for quote in ["Discussion: DCM.", "Experimental: DCE."]]
    arguments["source_conflicts"] = [{"description": "Discussion and procedure name different solvents", "excerpt_refs": refs}]
    event, block = prepare(workspace, arguments)
    assert set(refs) <= set(event.evidence_refs)
    assert any("source discrepancy" in note for note in block["limitations"])


def test_indexed_reuse_uses_saved_graphs_and_requires_same_publication(workspace):
    source, arguments = setup(workspace)
    lookup = workspace.store.append("call", {"operation": "get_precedents", "execution_status": "completed",
        "result": {"records": [{"reaction_id": "rxn1", "reference_id": PUBLICATION,
                                "reaction_smiles": "CCBr.N>>CCN"}]}})
    arguments.update(indexed_ref=lookup.artifact_ref, reaction_id="rxn1")
    for item in arguments["reactants"] + arguments["products"]:
        item.pop("smiles")
    event, block = prepare(workspace, arguments)
    assert block["structure_origin"] == "indexed_record"
    assert [p["smiles"] for p in block["reactants"]] == ["CCBr", "N"]
    assert lookup.artifact_ref in event.evidence_refs
    validate_answer_evidence(answer(source, block), workspace.store)
    different = workspace.capture_source(TEXT, url="https://example.org/paper", reference_id="REF1:" + "b" * 64)
    arguments["source_ref"] = different.artifact_ref
    failed = workspace.run("prepare_literature_reaction", arguments)
    assert "publication" in workspace.store.read_artifact(failed.artifact_ref)["error"]["message"].lower()


def test_answer_cannot_edit_prepared_graph_or_hide_source_limitations(workspace):
    source, arguments = setup(workspace)
    _, block = prepare(workspace, arguments)
    modified = deepcopy(block)
    modified["products"][0]["smiles"] = "CCC"
    with pytest.raises(ValueError, match="differs from"):
        validate_answer_evidence(answer(source, modified), workspace.store)
    modified = deepcopy(block)
    modified["source_provenance"]["acquisition"] = "http_fetch"
    with pytest.raises(ValueError, match="differs from"):
        validate_answer_evidence(answer(source, modified), workspace.store)


def test_source_file_import_preserves_bytes_and_bounded_cli_output(workspace, capsys):
    file = workspace.store.root / "browser.json"
    file.write_text(json.dumps({"document": {"text": TEXT}}), encoding="utf-8")
    result = main(["capture-source", str(workspace.store.root), "browser.json", "--url", "https://example.org/paper",
                   "--text-path", '["document","text"]', "--reference-id", PUBLICATION])
    assert result == 0
    output = capsys.readouterr().out
    assert TEXT not in output
    event_ref = json.loads(output)["event"]["artifact_ref"]
    record = workspace.store.read_artifact(event_ref)
    assert record["extraction"]["text"] == TEXT
    assert record["capture_file"]["text_path"] == ["document", "text"]
    with pytest.raises(ValueError, match="inside"):
        workspace.capture_source_file("../outside.txt", url="https://example.org/paper")


def test_caught_finalization_error_is_still_a_failed_activity(workspace):
    history = ActivityHistory()
    cursor = ScientificActivityCursor(workspace.store, after_sequence=0)
    with pytest.raises(ValueError):
        workspace.finalize_answer("answer-draft.json", {"answer_markdown": "Invalid claim", "claims": [{"id": "bad"}]})
    rows = cursor.drain(history)
    assert len(rows) == 1
    assert rows[0]["kind"] == "answer_validation"
    assert rows[0]["status"] == "failed"
    assert not (workspace.store.root / "answer-draft.json").exists()


def test_prepared_source_explicit_structures_and_svg_share_the_recorded_graph(workspace):
    text = "Compound 8: CCBr and ammonia N. Compound 9: CCN."
    source = workspace.capture_source(text, url="https://example.org/paper")
    passage = workspace.record_source_excerpt(source.artifact_ref, excerpt=text)
    _, arguments = setup(workspace)
    arguments.update(source_ref=source.artifact_ref, structure_evidence=text, structure_origin="source_explicit")
    for item in arguments["reactants"] + arguments["products"]:
        item["evidence_ref"] = passage.artifact_ref
        item.pop("reported_formula", None)
    _, block = prepare(workspace, arguments)
    payload = answer(source, block).model_dump()
    validate_answer_evidence(ScientificAnswer.model_validate(payload), workspace.store)
    view = present_conversation({"id": "a" * 32, "turns": [{"question": "Route?", "answer": payload}]})
    drawing = view["turns"][0]["structured_presentation"]["steps"][0]["literature_reactions"][0]
    assert drawing["drawing_status"] == "drawn"
    assert drawing["reaction_smiles"] == "CCBr.N>>CCN"
    assert drawing["source_provenance"]["acquisition"] == "agent_supplied_excerpt"
    assert drawing["products"][0]["evidence_artifact_url"].endswith(passage.artifact_ref)
    assert drawing["preparation_artifact_url"].endswith(block["preparation_ref"])


def test_forged_excerpt_and_invented_reported_conditions_are_rejected(workspace):
    source, arguments = setup(workspace)
    original = workspace.store.read_artifact(arguments["products"][0]["evidence_ref"])
    original["text"] = "Compound 9: invented passage"
    forged = workspace.store.append("literature_excerpt", original)
    arguments["products"][0]["evidence_ref"] = forged.artifact_ref
    failed = workspace.run("prepare_literature_reaction", arguments)
    assert "exactly" in workspace.store.read_artifact(failed.artifact_ref)["error"]["message"]
    _, arguments = setup(workspace)
    arguments["conditions"] = [{"text": "100 degrees for two hours", "basis": "reported",
                                 "source_ids": ["paper"], "limitations": []}]
    failed = workspace.run("prepare_literature_reaction", arguments)
    assert "quote" in workspace.store.read_artifact(failed.artifact_ref)["error"]["message"]
