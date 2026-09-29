"""Tool activity lifecycle and recovery from real CLI event shapes."""

import json

from chem_coworker.scientific_workspace.activity import ActivityHistory, activity_detail, recover_activity


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
