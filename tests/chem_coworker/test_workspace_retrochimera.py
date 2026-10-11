"""Optional supplier isolation, pinned dependencies, failure recording and reuse."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from chem_coworker.scientific_workspace import InvestigationStore, ScientificWorkspace
from chem_coworker.scientific_workspace.adapters import retrochimera as adapter
from chem_coworker.scientific_workspace.adapters.operations import ScientificOperations
from chem_coworker.scientific_workspace.adapters.retrochimera_config import load_config, prepare_config
from chem_coworker.scientific_workspace.core.baseline import artifact_identity, capture_baseline, verify_baseline

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def configuration(tmp_path):
    checkout = tmp_path / "external"
    source = checkout / "retrochimera"
    source.mkdir(parents=True)
    (source / "__init__.py").write_text("", "utf-8")
    model = checkout / "checkpoint"
    model.mkdir()
    (model / "models.json").write_text("{}", "utf-8")
    (model / "weights.bin").write_bytes(b"fake model")
    config = tmp_path / "supplier.json"
    config.write_text(json.dumps(prepare_config(checkout, Path(sys.executable), model)), "utf-8")
    return config


def test_missing_pins_and_new_dependencies_fail_closed(configuration):
    original = load_config(configuration)
    altered = {**original, "pinned_files": original["pinned_files"][:-1]}
    configuration.write_text(json.dumps(altered), "utf-8")
    with pytest.raises(ValueError, match="complete installation"):
        load_config(configuration)
    configuration.write_text(json.dumps(original), "utf-8")
    (Path(original["model_dir"]) / "new.bin").write_bytes(b"new model asset")
    with pytest.raises(ValueError, match="complete installation"):
        load_config(configuration)


def test_baseline_pins_checkpoint_contents(configuration):
    baseline = capture_baseline(ROOT, {"retrochimera_config": configuration})
    assert any(key.startswith("retrochimera_dependency_") for key in baseline["artifacts"])
    config = load_config(configuration)
    (Path(config["model_dir"]) / "weights.bin").write_bytes(b"changed weights")
    with pytest.raises(ValueError):
        verify_baseline(baseline)
    with pytest.raises(ValueError, match="checksum mismatch"):
        capture_baseline(ROOT, {"retrochimera_config": configuration})


@pytest.fixture
def fake_transport(tmp_path, monkeypatch):
    script = tmp_path / "worker.py"
    script.write_text('''import json, os, sys, time
from pathlib import Path
count = 0
for line in sys.stdin:
    directory = Path(json.loads(line))
    request = json.loads((directory / "request.json").read_text())
    if request.get("sleep"):
        time.sleep(60)
    if request.get("crash"):
        sys.exit(2)
    count += 1
    result = {"execution_status": "completed", "predictions": [],
              "execution": {"worker_pid": os.getpid(), "model_reused": count > 1}}
    temp = directory / "result.tmp"
    temp.write_text(json.dumps(result))
    temp.replace(directory / "result.json")
''', "utf-8")
    popen = subprocess.Popen

    def launch(args, **kwargs):
        args = [args[0], "-u", str(script), *args[3:]]
        return popen(args, **kwargs)

    monkeypatch.setattr(adapter.subprocess, "Popen", launch)


def test_persistent_worker_timeout_restart_and_cleanup(configuration, fake_transport, tmp_path):
    worker = adapter.RetroChimeraWorker()

    def request(index, **kwargs):
        directory = tmp_path / "diagnostics" / str(index)
        directory.mkdir(parents=True)
        (directory / "request.json").write_text(json.dumps(kwargs), "utf-8")
        return worker.request(configuration, directory, 0.3 if kwargs.get("sleep") else 10)

    try:
        first, second = request(1), request(2)
        assert first["execution_status"] == second["execution_status"] == "completed"
        assert not first["execution"]["model_reused"]
        assert second["execution"]["model_reused"]
        assert first["execution"]["worker_pid"] == second["execution"]["worker_pid"]
        process = worker.process
        timed_out = request(3, sleep=True)
        assert timed_out["execution_status"] == "timed_out"
        assert process.poll() is not None
        assert worker.process is None
        restarted = request(4)
        assert not restarted["execution"]["model_reused"]
        assert request(5, crash=True)["execution_status"] == "error"
    finally:
        worker.close()
    assert not worker.streams


def test_workspace_records_supplier_failure_without_success_claim(configuration, tmp_path, monkeypatch):
    store = InvestigationStore.create(tmp_path / "run", objective="Test supplier transport", baseline={
        "artifacts": {"retrochimera_config": artifact_identity(configuration)},
    })
    monkeypatch.setattr("chem_coworker.scientific_workspace.workspace.verify_baseline", lambda *_: None)
    monkeypatch.setattr(adapter.RetroChimeraWorker, "request", lambda *_: {
        "execution_status": "timed_out", "predictions": [],
        "error": {"type": "TimeoutExpired", "message": "test deadline"},
    })
    workspace = ScientificWorkspace(store.root)
    event = workspace.run("predict_retrochimera", {"target_smiles": "COC"})
    payload = store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "timed_out"
    assert payload["result"]["candidates"] == []
    assert payload["result"]["config_sha256"]
    assert "success" in " ".join(payload["result"]["limitations"])


def test_optional_supplier_missing_configuration_does_not_load_ml(tmp_path, monkeypatch):
    store = InvestigationStore.create(tmp_path / "run", objective="Optional supplier", baseline={"artifacts": {}})
    monkeypatch.setattr("chem_coworker.scientific_workspace.workspace.verify_baseline", lambda *_: None)
    workspace = ScientificWorkspace(store.root)
    event = workspace.run("predict_retrochimera", {"target_smiles": "COC"})
    payload = store.read_artifact(event.artifact_ref)
    assert payload["execution_status"] == "error"
    assert payload["error"]["type"] == "FileNotFoundError"
    assert ScientificOperations(store).definition("predict_retrochimera").required_artifacts == ("retrochimera_config",)


def test_replay_projection_keeps_model_identity_but_excludes_timing():
    definition = next(d for d in ScientificOperations.DEFINITIONS if d.name == "predict_retrochimera")
    first = {"candidates": [], "identity": {"config_sha256": "a"}, "execution": {"elapsed_seconds": 1}}
    second = {**first, "execution": {"elapsed_seconds": 7, "model_reused": True}}
    assert definition.replay_projection(first) == definition.replay_projection(second)
    assert definition.replay_projection(first) != definition.replay_projection({**second, "identity": {"config_sha256": "b"}})


def test_brief_view_exposes_candidates_and_reuse_without_changing_ranks():
    from chem_coworker.scientific_workspace.views.brief_summaries import summarize_call_brief

    candidates = [{"supplier_rank": i, "proposed_reaction_smiles": "CBr.CO>>COC",
                   "validation_status": "structures_parsed_reaction_unverified"} for i in range(2, 14)]
    summary = summarize_call_brief({"operation": "predict_retrochimera", "result": {
        "candidates": candidates, "execution": {"model_reused": True},
        "rejected": [{"supplier_rank": 1}],
    }})["result_summary"]
    assert summary["candidate_count"] == 12
    assert summary["candidates_omitted"] == 2
    assert summary["candidates"][0]["supplier_rank"] == 2
    assert summary["rejected_count"] == 1
    assert summary["execution"]["model_reused"] is True
