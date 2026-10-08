"""Tool activity lifecycle and recovery from real CLI event shapes."""

import json

import pytest

from chem_coworker.scientific_workspace.runtime.activity import (
    ActivityHistory,
    activity_detail,
    recover_activity,
)


@pytest.mark.parametrize("receipt", [
    {"exit_code": 1, "stdout": "", "stderr": "ValueError: invalid request"},
    {"returncode": 2, "stdout": "", "stderr": "ValueError: invalid request"},
    {"error": "Command failed", "out": "", "err": "ValueError: invalid request"},
])
@pytest.mark.parametrize("structured", [False, True])
def test_nested_process_failure_is_visible_and_recoverable(tmp_path, receipt, structured):
    result = {"structured_content": receipt} if structured else {
        "content": [{"type": "text", "text": json.dumps(receipt)}],
    }
    event = {"type": "item.completed", "item": {
        "id": "node", "type": "mcp_tool_call", "server": "node_repl", "tool": "js",
        "status": "completed", "result": result,
    }}
    row = ActivityHistory().observe(event, "now")
    assert row["status"] == "failed"
    assert "ValueError: invalid request" in row["failure_detail"]
    path = tmp_path / "runtime.jsonl"
    path.write_text(json.dumps(event), encoding="utf-8")
    before = path.read_bytes()
    assert recover_activity([path], [])[0]["status"] == "failed"
    assert path.read_bytes() == before


@pytest.mark.parametrize("text", [
    'Example traceback: ValueError: invalid request',
    '{"error": "unresolved chemistry", "status": "review"}',
    '{"exit_code": 0, "stdout": "ok", "stderr": "warning"}',
    '',
])
def test_node_output_is_not_guessed_to_be_a_process_failure(text):
    row = ActivityHistory().observe({"type": "item.completed", "item": {
        "id": "node", "type": "mcp_tool_call", "server": "node_repl", "tool": "js",
        "status": "completed", "result": {"content": [{"type": "text", "text": text}]},
    }}, "now")
    assert row["status"] == "completed"
    assert not row["failure_detail"]


def test_nonzero_exit_overrides_outer_completed_event():
    row = ActivityHistory().observe({"type": "item.completed", "item": {
        "id": "shell", "type": "command_execution", "status": "completed",
        "exit_code": 1, "aggregated_output": "ValueError: invalid request",
    }}, "now")
    assert row["status"] == "failed"
    assert row["failure_detail"] == "ValueError: invalid request"


@pytest.mark.parametrize("prefix", [
    "pydantic_core._pydantic_core.ValidationError: 10 validation errors for ScientificAnswer\n",
    "pydantic_core.ValidationError: 10 validation errors for ScientificAnswer\n" + "x" * 5000 + "\n",
])
def test_validation_failure_preserves_field_and_message_without_documentation_footer(prefix):
    output = (prefix + "steps.0.conditions.0.text\n"
              "  Field required [type=missing, input_value={}, input_type=dict]\n"
              "    For further information visit https://errors.pydantic.dev/2.11/v/missing\n")
    history = ActivityHistory()
    row = history.observe({"type": "item.completed", "item": {
        "id": "validation", "type": "command_execution", "exit_code": 1,
        "aggregated_output": output,
    }}, "finish")
    assert "steps.0.conditions.0.text" in row["failure_detail"]
    assert "Field required" in row["failure_detail"]
    assert "https://errors.pydantic.dev" not in row["failure_detail"]
    assert len(row["failure_detail"]) <= 500


def test_regular_traceback_keeps_final_exception_and_success_hides_output():
    history = ActivityHistory()
    for code, expected in [(1, "ValueError: invalid input"), (0, "")]:
        row = history.observe({"type": "item.completed", "item": {
            "id": str(code), "type": "command_execution", "exit_code": code,
            "aggregated_output": "Traceback (most recent call last):\nValueError: invalid input\n",
        }}, "finish")
        assert row["failure_detail"] == expected


def test_command_lifecycle_keeps_action_and_bounded_failure_in_one_row():
    history = ActivityHistory()
    history.observe({"type": "item.started", "item": {
        "id": "cmd", "type": "command_execution", "status": "in_progress",
        "command": '"C:\\Program Files\\PowerShell\\7\\pwsh.exe" -Command "Get-Content docs/old.md"',
    }}, "start")
    history.observe({"type": "item.completed", "item": {
        "id": "cmd", "type": "command_execution", "exit_code": 1,
        "aggregated_output": "x" * 9000 + "\n\x1b[31mCannot find path docs/old.md because it does not exist.\x1b[0m",
    }}, "finish")
    assert len(history.rows) == 1
    row = history.rows[0]
    assert row["title"] == "Read old.md"
    assert row["at"] == "start" and row["updated_at"] == "finish"
    assert row["exit_code"] == 1
    assert "Cannot find path" in row["failure_detail"]
    assert len(row["failure_detail"]) <= 500
    assert "aggregated_output" not in row
    assert "\x1b" not in row["failure_detail"]


def test_query_arriving_at_completion_replaces_pending_search():
    history = ActivityHistory()
    history.observe({"type": "item.started", "item": {
        "id": "search", "type": "web_search", "query": "", "action": {"type": "other"},
    }}, "start")
    assert "waiting for query details" in history.rows[0]["title"]
    history.observe({"type": "item.updated", "item": {
        "id": "search", "action": {"type": "search", "queries": ["exact target synthesis", "patent Example 2"]},
    }}, "update")
    history.observe({"type": "item.completed", "item": {"id": "search", "type": "web_search", "query": ""}}, "end")
    assert len(history.rows) == 1
    assert history.rows[0]["detail"] == "exact target synthesis; patent Example 2"
    assert history.rows[0]["status"] == "completed"
    assert "waiting" not in history.rows[0]["title"]


def test_completed_web_action_without_metadata_is_not_still_waiting_or_verified():
    history = ActivityHistory()
    for phase in ("started", "completed"):
        history.observe({"type": "item." + phase, "item": {
            "id": "web", "type": "web_search", "query": "", "action": {"type": "other"},
        }}, phase)
    row = history.rows[0]
    assert len(history.rows) == 1
    assert row["status"] == "completed" and "waiting" not in row["title"]
    assert "unavailable" in row["title"]
    assert row["result_availability"] == "not_reported" and row["source_url"] is None


def test_web_error_and_known_source_are_preserved_even_on_completion():
    history = ActivityHistory()
    history.observe({"type": "item.completed", "item": {
        "id": "web", "type": "web_search", "action": {
            "type": "screenshot", "url": "https://example.org/paper.pdf", "pageno": 0,
        }, "error": {"message": "No usable image"},
    }}, "end")
    row = history.rows[0]
    assert row["status"] == "failed" and row["failure_detail"] == "No usable image"
    assert row["source_url"] == "https://example.org/paper.pdf"
    assert row["web_action_type"] == "screenshot" and row["result_availability"] == "failed"
    assert "page 1" in row["detail"]


def test_retry_ids_are_distinct_and_success_output_and_reasoning_are_excluded():
    history = ActivityHistory()
    event = {"type": "item.completed", "item": {
        "id": "item_1", "type": "command_execution", "command": "python analyze.py",
        "exit_code": 0, "aggregated_output": "Large private result",
    }}
    history.observe(event, "first", scope="0")
    history.observe(event, "second", scope="1")
    history.observe({"type": "item.completed", "item": {"id": "reasoning", "type": "reasoning", "text": "Hidden"}}, "third")
    assert len(history.rows) == 2
    assert history.rows[0]["activity_id"] != history.rows[1]["activity_id"]
    assert history.rows[0]["title"] == "Run Python script: analyze.py"
    assert "Large private result" not in json.dumps(history.rows)
    assert all(not row["failure_detail"] for row in history.rows)


def test_web_page_and_scientific_tool_details_use_observed_inputs():
    assert activity_detail({"type": "web_search", "action": {
        "type": "open_page", "url": "https://example.org/paper",
    }}) == "Open page: https://example.org/paper"
    assert activity_detail({"type": "web_search", "action": {
        "type": "find_in_page", "url": "https://example.org/paper", "pattern": "yield",
    }}) == "Find yield in https://example.org/paper"
    detail = activity_detail({"type": "mcp_tool_call", "server": "chemistry", "tool": "disconnect_target",
                              "arguments": '{"target_smiles":"CCN","unrelated":"do not display"}'})
    assert "disconnect_target" in detail and "target_smiles: CCN" in detail
    assert "unrelated" not in detail


@pytest.mark.parametrize("encoded", [False, True])
def test_node_repl_title_code_output_and_duration_survive_lifecycle_and_recovery(tmp_path, encoded):
    arguments = {"title": "Inspect scientific workspace", "code": 'const x = "<tag>&amp;";\nnodeRepl.write(x);'}
    events = [
        {"type": "item.started", "item": {"id": "node", "type": "mcp_tool_call",
         "server": "node_repl", "tool": "js", "arguments": json.dumps(arguments) if encoded else arguments}},
        {"type": "item.completed", "item": {"id": "node", "type": "mcp_tool_call", "result": {
            "content": [{"type": "text", "text": "start\n" + "x" * 9000 + "\nValueError: nested process failed"},
                        {"type": "image", "data": "image bytes excluded"}],
            "_meta": {"codex/nodeReplExecutionDurationMs": 54335}}}},
    ]
    history = ActivityHistory()
    for event in events:
        history.observe(event, "now")
    row = history.rows[0]
    assert len(history.rows) == 1
    assert row["action_title"] == row["title"] == arguments["title"]
    assert row["code_preview"] == arguments["code"]
    assert row["detail"] == "node_repl / js"
    assert row["duration_ms"] == 54335
    assert len(row["output_preview"]) <= 2400
    assert row["output_preview"].startswith("start\n")
    assert row["output_preview"].endswith("ValueError: nested process failed")
    assert "truncated" in row["output_preview"]
    assert "image bytes" not in json.dumps(row)
    # Output may describe a nested failure; do not guess the outer tool's status.
    assert row["status"] == "completed" and not row["failure_detail"]
    path = tmp_path / "runtime.jsonl"
    path.write_text("\n".join(json.dumps(event) for event in events), encoding="utf-8")
    original = path.read_bytes()
    recovered = recover_activity([path], [])
    assert recovered[0]["action_title"] == arguments["title"]
    assert recovered[0]["output_preview"] == row["output_preview"]
    assert path.read_bytes() == original
    row = history.observe({"type": "item.completed", "item": {"id": "node"}}, "later")
    assert row["duration_ms"] == 54335 and row["output_preview"] == recovered[0]["output_preview"]


@pytest.mark.parametrize("arguments", [None, "{invalid", "[]", {"title": 123}, {"title": "  "}])
def test_mcp_without_valid_title_keeps_tool_fallback(arguments):
    row = ActivityHistory().observe({"type": "item.completed", "item": {
        "id": "tool", "type": "mcp_tool_call", "server": "node_repl", "tool": "js",
        "arguments": arguments,
    }}, "now")
    assert row["action_title"] == ""
    assert row["detail"] == "node_repl / js"


def test_mcp_error_result_is_failure_and_successful_other_tools_do_not_expose_output():
    for failed in (False, True):
        row = ActivityHistory().observe({"type": "item.completed", "item": {
            "id": "tool", "type": "mcp_tool_call", "server": "chemistry", "tool": "analyze",
            "result": {"isError": failed, "content": [{"type": "text", "text": "Diagnostic text"}]},
        }}, "now")
        assert row["status"] == ("failed" if failed else "completed")
        assert row["failure_detail"] == ("Diagnostic text" if failed else "")
        assert "output_preview" not in row


def test_recovery_restores_dropped_history_without_rewriting_logs(tmp_path):
    path = tmp_path / "runtime.jsonl"
    events = [{"type": "item.completed", "item": {
        "id": f"item_{i}", "type": "command_execution", "status": "completed",
        "command": f"python analysis_{i}.py", "exit_code": 0,
    }} for i in range(40)]
    path.write_text("\n".join(json.dumps(event) for event in events) + '\n{"partial":', encoding="utf-8")
    before = path.read_bytes()
    saved = [{"kind": "command_execution", "status": "completed", "at": str(i)} for i in range(10, 40)]
    rows = recover_activity([path, tmp_path / "missing.jsonl"], saved)
    assert len(rows) == 40
    assert rows[0]["title"] == "Run Python script: analysis_0.py"
    assert rows[0]["at"] is None
    assert rows[10]["at"] == "10" and rows[-1]["at"] == "39"
    assert path.read_bytes() == before
    assert recover_activity([tmp_path / "missing.jsonl"], saved) == saved


def test_recovery_keeps_scientific_outcomes_between_runtime_actions(tmp_path):
    path = tmp_path / "runtime.jsonl"
    events = [{"type": "item.completed", "item": {
        "id": str(i), "type": "mcp_tool_call", "status": "completed", "server": "node_repl",
        "tool": "js", "arguments": {"title": f"Action {i}"},
    }} for i in (1, 3)]
    path.write_text("\n".join(json.dumps(event) for event in events), encoding="utf-8")
    scientific = {"kind": "scientific_call", "activity_id": "scientific:4", "at": "2",
                  "status": "failed", "failure_detail": "Invalid input", "event_sequence": 4}
    saved = [{"kind": "mcp_tool_call", "status": "completed", "at": "1"}, scientific,
             {"kind": "mcp_tool_call", "status": "completed", "at": "3"}]
    rows = recover_activity([path], saved)
    assert [row["at"] for row in rows] == ["1", "2", "3"]
    assert rows[1] == scientific
    assert rows[0]["action_title"] == "Action 1" and rows[2]["action_title"] == "Action 3"


def test_only_public_commentary_is_shown_and_recovered(tmp_path):
    history = ActivityHistory()
    messages = [
        {"id": "comment", "type": "agent_message", "text": "I found a reported route; I’m checking its yield."},
        {"id": "private", "type": "reasoning", "text": "Private reasoning"},
        {"id": "phase", "type": "agent_message", "phase": "analysis", "text": "Private analysis"},
        {"id": "final", "type": "agent_message", "phase": "final_answer", "text": "Final answer"},
        {"id": "json", "type": "agent_message", "text": '{"schema_version":"scientific_answer.v2"}'},
        {"id": "fenced", "type": "agent_message", "text": '```json\n{"answer_markdown":"Final"}\n```'},
    ]
    events = [{"type": "item.completed", "item": item} for item in messages]
    for event in events:
        history.observe(event, "now")
    assert len(history.rows) == 1
    assert history.rows[0]["kind"] == "agent_update"
    assert history.rows[0]["detail"] == messages[0]["text"]
    path = tmp_path / "runtime.jsonl"
    path.write_text("\n".join(json.dumps(event) for event in events), "utf-8")
    recovered = recover_activity([path], [])
    assert len(recovered) == 1 and recovered[0]["detail"] == messages[0]["text"]
