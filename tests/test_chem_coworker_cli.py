import json
from pathlib import Path

import pytest

from chem_coworker._cli import app as cli_app
from chem_coworker._cli import commands as cli_commands
from chem_coworker._cli.commands import (
    COMMAND_EXIT,
    COMMAND_HANDLED,
    COMMAND_UNHANDLED,
    ReplSession,
    build_default_command_registry,
)
from chem_coworker.response import ChemResponse


class _FakeCoworker:
    def __init__(self) -> None:
        self.verbose = False
        self.plan_callback = None
        self.compact_calls = 0
        self.run_queries = []

    def set_verbose(self, verbose: bool) -> None:
        self.verbose = verbose

    def set_plan_callback(self, callback) -> None:
        self.plan_callback = callback

    def compact_history(self, history):
        self.compact_calls += 1
        return [{"role": "assistant", "content": "summary"}]

    def run(self, query: str) -> ChemResponse:
        self.run_queries.append(query)
        if query.startswith("FAIL:"):
            raise RuntimeError(query)
        return ChemResponse(
            query=query,
            task_type="general",
            answer=f"Answer: {query}",
            tools_called=[],
            tool_results={},
            model="fake-model",
            provider="openai",
            elapsed_s=0.1,
            llm_calls=1,
            streamed=False,
        )


class _FakeToolRegistry:
    def describe_tools(self) -> str:
        return "fake tools"


class _FakeUI:
    def __init__(self) -> None:
        self.settings_calls = []

    def show_plan_and_confirm(self, plan):
        return plan

    def print_settings(self, model: str, provider: str, verbose: bool, plan_mode: bool) -> None:
        self.settings_calls.append((model, provider, verbose, plan_mode))


def _make_session() -> ReplSession:
    coworker = _FakeCoworker()
    ui = _FakeUI()
    return ReplSession(
        coworker=coworker,
        model="o4-mini",
        provider="openai",
        verbose=False,
        plan_mode=False,
        history=[{"role": "user", "content": "x"}, {"role": "assistant", "content": "y"}],
        tool_registry=_FakeToolRegistry(),
        ui=ui,
        create_coworker=lambda m, p, v, pm: coworker,
    )


def test_command_registry_toggles_verbose_and_plan(capsys) -> None:
    reg = build_default_command_registry()
    session = _make_session()

    assert reg.dispatch("/verbose", session) == COMMAND_HANDLED
    assert session.verbose is True
    assert session.coworker.verbose is True

    assert reg.dispatch("/plan", session) == COMMAND_HANDLED
    assert session.plan_mode is True
    assert session.coworker.plan_callback is not None

    assert reg.dispatch("/history", session) == COMMAND_HANDLED
    out = capsys.readouterr().out
    assert "History:" in out


def test_command_registry_compact_and_exit() -> None:
    reg = build_default_command_registry()
    session = _make_session()

    assert reg.dispatch("/compact", session) == COMMAND_HANDLED
    assert session.coworker.compact_calls == 1
    assert session.history == [{"role": "assistant", "content": "summary"}]

    assert reg.dispatch("quit", session) == COMMAND_EXIT
    assert reg.dispatch("/not-a-command", session) == COMMAND_UNHANDLED


def test_command_registry_session_save_list_load(monkeypatch, tmp_path, capsys) -> None:
    reg = build_default_command_registry()
    session = _make_session()
    sessions_dir = tmp_path / "sessions"
    monkeypatch.setattr("chem_coworker._cli.session_store.SESSIONS_DIR", sessions_dir)

    assert reg.dispatch("/session save demo", session) == COMMAND_HANDLED
    assert (sessions_dir / "demo.json").exists()

    out = capsys.readouterr().out
    assert "Session saved" in out

    assert reg.dispatch("/session list", session) == COMMAND_HANDLED
    out = capsys.readouterr().out
    assert "demo" in out

    session.history = []
    assert reg.dispatch("/session load demo", session) == COMMAND_HANDLED
    assert len(session.history) == 2
    out = capsys.readouterr().out
    assert "Session loaded" in out


def test_parser_help_includes_ask_and_batch() -> None:
    parser = cli_app._build_parser()
    help_text = parser.format_help()
    assert "ask" in help_text
    assert "batch" in help_text


def test_main_ask_json_outputs_chemresponse(monkeypatch, capsys) -> None:
    fake = _FakeCoworker()
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setattr(cli_app, "_init_coworker", lambda *args, **kwargs: fake)
    monkeypatch.setattr(cli_app, "_resolve_model_provider", lambda args: ("o4-mini", "openai"))

    cli_app.main(["ask", "hello world", "--json"])
    out = capsys.readouterr().out
    payload = json.loads(out)
    assert payload["query"] == "hello world"
    assert payload["answer"] == "Answer: hello world"
    assert payload["model"] == "fake-model"


def test_main_ask_output_format_json_outputs_chemresponse(monkeypatch, capsys) -> None:
    fake = _FakeCoworker()
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setattr(cli_app, "_init_coworker", lambda *args, **kwargs: fake)
    monkeypatch.setattr(cli_app, "_resolve_model_provider", lambda args: ("o4-mini", "openai"))

    cli_app.main(["ask", "hello", "--output-format", "json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["query"] == "hello"


def test_main_batch_jsonl_from_lines(monkeypatch, capsys, tmp_path) -> None:
    fake = _FakeCoworker()
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setattr(cli_app, "_init_coworker", lambda *args, **kwargs: fake)
    monkeypatch.setattr(cli_app, "_resolve_model_provider", lambda args: ("o4-mini", "openai"))

    batch_file = tmp_path / "queries.txt"
    batch_file.write_text("# comment\nfirst query\n\nsecond query\n", encoding="utf-8")

    cli_app.main(["batch", str(batch_file), "--json"])
    out_lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.strip()]
    assert len(out_lines) == 2

    rec1 = json.loads(out_lines[0])
    rec2 = json.loads(out_lines[1])
    assert rec1["index"] == 1
    assert rec1["query"] == "first query"
    assert rec2["index"] == 2
    assert rec2["query"] == "second query"


def test_main_batch_jsonl_failure_exits_nonzero_and_writes_summary(monkeypatch, capsys, tmp_path) -> None:
    fake = _FakeCoworker()
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setattr(cli_app, "_init_coworker", lambda *args, **kwargs: fake)
    monkeypatch.setattr(cli_app, "_resolve_model_provider", lambda args: ("o4-mini", "openai"))

    batch_file = tmp_path / "queries.txt"
    batch_file.write_text("ok one\nFAIL:boom\nok two\n", encoding="utf-8")

    with pytest.raises(SystemExit) as excinfo:
        cli_app.main(["batch", str(batch_file), "--output-format", "jsonl"])
    assert excinfo.value.code == 1

    stdout_lines = [ln for ln in capsys.readouterr().out.splitlines() if ln.strip()]
    # 3 JSONL records should be emitted despite one failure (continue-on-error default)
    assert len(stdout_lines) == 3
    recs = [json.loads(ln) for ln in stdout_lines]
    assert recs[1]["success"] is False
    # Summary is emitted to stderr as JSON to avoid corrupting JSONL stdout
    summary_json = capsys.readouterr().err
    assert summary_json == ""  # capsys has already been consumed above


def test_main_batch_jsonl_failure_summary_in_stderr(monkeypatch, capsys, tmp_path) -> None:
    fake = _FakeCoworker()
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setattr(cli_app, "_init_coworker", lambda *args, **kwargs: fake)
    monkeypatch.setattr(cli_app, "_resolve_model_provider", lambda args: ("o4-mini", "openai"))

    batch_file = tmp_path / "queries.txt"
    batch_file.write_text("FAIL:boom\n", encoding="utf-8")

    with pytest.raises(SystemExit):
        cli_app.main(["batch", str(batch_file), "--output-format", "jsonl", "--fail-fast"])
    captured = capsys.readouterr()
    assert captured.out.strip()
    summary = json.loads(captured.err.strip())
    assert summary["type"] == "batch_summary"
    assert summary["failed"] == 1


def test_main_batch_output_file_requires_jsonl(monkeypatch, tmp_path) -> None:
    fake = _FakeCoworker()
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    monkeypatch.setattr(cli_app, "_init_coworker", lambda *args, **kwargs: fake)
    monkeypatch.setattr(cli_app, "_resolve_model_provider", lambda args: ("o4-mini", "openai"))

    batch_file = tmp_path / "queries.txt"
    batch_file.write_text("ok one\n", encoding="utf-8")
    out_file = tmp_path / "out.txt"

    with pytest.raises(SystemExit) as excinfo:
        cli_app.main(["batch", str(batch_file), "--output", str(out_file)])
    assert excinfo.value.code == 2


def test_main_intake_json_dry_run_outputs_payload(monkeypatch, capsys) -> None:
    class _FakeExtractor:
        def __init__(self, provider=None, model=None, verbose=False):  # noqa: ARG002
            pass

        def intake(self, **kwargs):
            assert kwargs["dry_run"] is True
            assert kwargs["mismatch_policy"] == "warn"
            return {
                "success": True,
                "source": "demo",
                "note_type": "reactions",
                "reaction_types": ["suzuki_miyaura"],
                "reaction_type_hint_raw": "",
                "reaction_type_hint_canonical": None,
                "reaction_types_detected_raw": ["suzuki"],
                "reaction_types_detected_canonical": ["suzuki_miyaura"],
                "reaction_types_unknown": [],
                "reaction_type_suggestions": {},
                "mismatch": False,
                "mismatch_policy": "warn",
                "unknown_reaction_policy": "general",
                "dry_run": True,
                "write_performed": False,
                "quarantined": False,
                "quarantine_file": "",
                "notes_files": ["knowledge_base/notes/reactions/suzuki_miyaura.md"],
                "extracted_notes": "notes",
                "char_count": 5,
                "warnings": [],
                "warning_details": [],
            }

    monkeypatch.setattr("chem_coworker.extractor.NotesExtractor", _FakeExtractor)
    monkeypatch.setattr(cli_app, "load_config", lambda: {"name": "o4-mini", "provider": "openai"})

    cli_app.main(["intake", "demo.txt", "--dry-run", "--output-format", "json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["success"] is True
    assert payload["dry_run"] is True
    assert payload["reaction_types"] == ["suzuki_miyaura"]


def test_main_intake_list_reaction_types_json(monkeypatch, capsys) -> None:
    monkeypatch.setattr(
        "chemtools.taxonomy.reaction_catalog.list_reaction_type_ids",
        lambda: ["suzuki_miyaura", "heck"],
    )
    cli_app.main(["intake", "--list-reaction-types", "--output-format", "json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["count"] == 2
    assert payload["reaction_type_ids"][0] == "suzuki_miyaura"


def test_load_batch_queries_supports_jsonl_string_and_object(tmp_path) -> None:
    path = tmp_path / "queries.jsonl"
    path.write_text(
        json.dumps("a") + "\n" + json.dumps({"query": "b"}) + "\n",
        encoding="utf-8",
    )
    assert cli_app._load_batch_queries(str(path), "auto") == ["a", "b"]
