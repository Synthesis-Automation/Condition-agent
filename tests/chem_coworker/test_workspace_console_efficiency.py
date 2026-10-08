"""Bounded console views retain literature leads and reuse complete saved evidence."""

import json

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.__main__ import main
from chem_coworker.scientific_workspace.adapters import literature
from chem_coworker.scientific_workspace.core.store import canonical_bytes


def test_fragment_brief_retains_publication_and_graph_without_raw_columns(tmp_path, capsys):
    store = InvestigationStore.create(tmp_path / "study", objective="Inspect source", baseline={})
    raw = {"authors": "Ariyasu, Kosuke", "citation": "American Journal of Organic Chemistry (2016), 6(3), 93-101",
           "irrelevant_large_column": "x" * 200000}
    record = {"reaction_smiles": "CCBr.N>>CCN", "yield_pct": 90,
              "source": {"raw_fields": {"source_observation": {"raw_fields": raw}}}}
    hit = {"reaction_id": "rxn1", "reference_id": "REF1:" + "a" * 64, "record": record,
           "product_smiles": "CCN", "warnings": ["SOURCE_MAPPING_UNRESOLVED"]}
    event = store.append("call", {"operation": "search_fragment_precedents", "execution_status": "completed",
                                  "result": {"hits": [hit] * 6, "search_status": "completed"}})
    workspace = ScientificWorkspace(store.root)
    before = store.read_artifact(event.artifact_ref)
    count = len(store.events())
    summary = workspace.call_summary(event.artifact_ref)
    first = summary["result_summary"]["hits"][0]
    assert first["publication"]["authors"] == "Ariyasu, Kosuke"
    assert "2016" in first["publication"]["citation"]
    assert first["reaction_smiles"] == "CCBr.N>>CCN"
    assert first["warnings"] == ["SOURCE_MAPPING_UNRESOLVED"]
    assert summary["result_summary"]["hit_page"]["next_offset"] == 3
    assert len(canonical_bytes(summary)) < 6000
    assert main(["call-summary", str(store.root), event.artifact_ref]) == 0
    assert "irrelevant_large_column" not in capsys.readouterr().out
    assert main(["show", str(store.root), event.artifact_ref]) == 0
    bounded = json.loads(capsys.readouterr().out)
    assert bounded["inspection_schema_version"] == "scientific_artifact_inspection.v1"
    assert main(["show", str(store.root), event.artifact_ref, "--full"]) == 0
    assert json.loads(capsys.readouterr().out) == before
    assert store.read_artifact(event.artifact_ref) == before
    assert len(store.events()) == count  # Reads do not recompute or fabricate events.


def test_oversized_summary_cannot_dump_large_structures_or_hide_omissions(tmp_path):
    store = InvestigationStore.create(tmp_path / "study", objective="Large route", baseline={})
    steps = [{"external_step_id": f"s{i}", "assessment": {"status": "unsupported", "warnings": ["CONFLICT"],
             "canonical_target_smiles": "C" * 4000, "canonical_precursor_smiles": "C" * 4000}}
             for i in range(20)]
    event = store.append("call", {"operation": "assess_route_proposal", "execution_status": "completed",
                                  "result": {"assessment": {"status": "partially_supported", "step_assessments": steps}}})
    summary = ScientificWorkspace(store.root).call_summary(event)
    assert len(canonical_bytes(summary)) <= 16000
    if summary["inspection"].get("summary_omitted"):
        assert "Inspect" in summary["inspection"]["hint"]
    else:
        assert summary["result_summary"]["assessment"]["status"] == "partially_supported"


def test_batch_budget_is_shared_and_omitted_references_remain_inspectable(tmp_path):
    store = InvestigationStore.create(tmp_path / "study", objective="Batch inspection", baseline={})
    events = [store.append("call", {"operation": "analyze_molecule", "execution_status": "completed",
        "result": {"valid": True, "canonical_smiles": "C" * 2000, "warnings": ["x" * 180] * 30,
                   "motifs": [{"motif_id": str(i), "chemist_label": "z" * 180}] * 30}}) for i in range(30)]
    workspace = ScientificWorkspace(store.root)
    before = len(store.events())
    summary = workspace.batch_summary([e.artifact_ref for e in events])
    assert len(canonical_bytes(summary)) <= 16000
    assert any(item.get("summary_omitted") for item in summary["items"])
    assert {item.get("artifact_ref", item.get("event", {}).get("artifact_ref")) for item in summary["items"]} == {e.artifact_ref for e in events}
    assert len(store.events()) == before
    help_entries = workspace.help(["capture_source_file", "assess_proposed_recipe"])
    assert "never a DOI" in help_entries[0]["description"]
    assert "operating_conditions" in help_entries[1]["signature"]
    entries = workspace.help(["run", "run_summary", "call_summary"])
    assert [entry["name"] for entry in entries] == ["run", "run_summary", "call_summary"]
    assert all(entry["signature"] and entry["example"] for entry in entries)
    assert len(store.events()) == before


def test_network_denial_is_not_repeated_without_explicit_transport_change(tmp_path, monkeypatch):
    store = InvestigationStore.create(tmp_path / "study", objective="Capture paper", baseline={})
    calls = []

    def denied(url, **kwargs):
        calls.append(url)
        raise PermissionError(13, "Fixture: network permission denied")

    monkeypatch.setattr(literature, "_request", denied)
    first = literature.fetch_source(store, "https://example.org/paper")
    second = literature.fetch_source(store, "https://example.org/paper.pdf")
    assert len(calls) == 1
    saved = store.read_artifact(second.artifact_ref)
    assert saved["network_request_skipped"] is True
    assert saved["blocked_by_ref"] == first.artifact_ref
    assert saved["source_url"].endswith("paper.pdf")
    assert saved["recovery"]["source_url"].endswith("paper.pdf")
    literature.fetch_source(store, "https://example.org/paper.pdf", retry_network=True)
    assert len(calls) == 2
