"""Readable tool activity projected from runtime events, never model reasoning."""

from __future__ import annotations

from html import unescape
import json
from pathlib import Path
import re
from typing import Any, Mapping


ACTIVITY_VERSION = 3
_EVENTS = {"item.started", "item.updated", "item.completed", "item.failed"}
_KINDS = {"command_execution", "web_search", "mcp_tool_call", "file_change", "todo_list", "agent_message"}


def _text(value: Any, limit: int = 600) -> str:
    if not isinstance(value, str):
        return ""
    value = re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", unescape(value))
    value = " ".join(value.split())
    return value if len(value) <= limit else value[:limit - 1] + "…"


def activity_detail(item: Mapping[str, Any]) -> str:
    """Show actual commands, search actions, tool arguments or changed files."""
    kind = item.get("type")
    if kind == "agent_message":
        return _text(item.get("text"), 1200)
    if kind == "command_execution":
        return _text(item.get("command"))
    if kind == "web_search":
        action = item.get("action")
        action = action if isinstance(action, Mapping) else {}
        queries = action.get("queries")
        if isinstance(queries, list) and queries:
            return _text("; ".join(query for query in queries if isinstance(query, str)))
        if action.get("type") == "open_page":
            return _text(f"Open page: {action.get('url', '')}")
        if action.get("type") == "find_in_page":
            return _text(f"Find {action.get('pattern', '')} in {action.get('url', '')}")
        return _text(item.get("query"))
    if kind == "mcp_tool_call":
        name = " / ".join(str(item[key]) for key in ("server", "tool") if item.get(key))
        arguments = item.get("arguments")
        if isinstance(arguments, str):
            try:
                arguments = json.loads(arguments)
            except ValueError:
                arguments = {}
        if isinstance(arguments, Mapping):
            selected = [f"{key}: {_text(arguments[key], 160)}" for key in (
                "target_smiles", "reaction_smiles", "smiles", "query", "url", "operation", "source_ref",
            ) if isinstance(arguments.get(key), str)]
            if selected:
                name += " — " + "; ".join(selected)
        return _text(name)
    if kind == "file_change":
        return _text("; ".join(
            f"{change.get('kind', 'update')}: {change['path']}"
            for change in (item.get("changes") or []) if isinstance(change, dict) and change.get("path")
        ))
    if kind == "todo_list":
        return _text("; ".join(
            entry["text"] for entry in (item.get("items") or [])
            if isinstance(entry, dict) and isinstance(entry.get("text"), str) and not entry.get("completed")
        )) or "Plan updated"
    return ""


def _command_title(command: str) -> str:
    # Strip the shell launcher, while leaving the actual command in the detail.
    body = re.split(r"\s+(?:-Command|-c|-lc)\s+", command, maxsplit=1, flags=re.I)[-1].strip('"')
    operation = re.search(r"\.run\(\s*['\"]([a-z_]+)['\"]", body)
    if operation:
        return "Run scientific operation: " + operation[1]
    script = re.search(r"(?:^|[\s'\"\\/])([^\s'\"\\/]+\.py)\b", body)
    if script and re.search(r"python|py\.exe", body, re.I):
        return "Run Python script: " + script[1]
    read = re.search(r"(?:Get-Content|\bcat)\s+(?:(?:-LiteralPath|-Path)\s+)?['\"]?([^\s'\";]+)", body, re.I)
    if read:
        return "Read " + re.split(r"[\\/]", read[1])[-1]
    if re.match(r"(?:rg|Select-String)\b", body, re.I):
        return "Search local files: " + _text(body, 120)
    if re.match(r"(?:Get-Location|Get-ChildItem|pwd|ls)\b", body, re.I):
        return "Inspect workspace files and paths"
    return "Run command: " + _text(body, 140)


def _title(item: Mapping[str, Any], detail: str) -> str:
    kind = item.get("type")
    if kind == "command_execution":
        return _command_title(_text(item.get("command"), 4000)) if detail else "Run command"
    if kind == "web_search":
        action = item.get("action") or {}
        action_type = action.get("type") if isinstance(action, Mapping) else None
        label = {"open_page": "Read web page", "find_in_page": "Find text on web page"}.get(action_type, "Search the web")
        if action_type in {None, "other"} and detail.startswith(("https://", "http://")):
            label = "Web request"
        return label + (": " + _text(detail, 150) if detail else " — waiting for query details")
    if kind == "mcp_tool_call":
        return "Call " + _text(item.get("tool") or "scientific tool", 150)
    return {"file_change": "Update investigation files", "todo_list": "Update investigation plan"}.get(kind, "Agent activity")


def _failure(item: Mapping[str, Any]) -> str:
    code = item.get("exit_code")
    if item.get("status") != "failed" and not (type(code) is int and code != 0):
        return ""
    error = item.get("error")
    if isinstance(error, Mapping):
        error = error.get("message")
    if _text(error):
        return _text(error, 500)
    # Only failed actions expose a bounded diagnostic tail, never successful output.
    output = item.get("aggregated_output")
    if isinstance(output, str) and output.strip():
        lines = [_text(line, 500).lstrip("| ") for line in output[-4000:].splitlines()]
        lines = [line for line in lines if line and not re.fullmatch(r"[~^\-\s]+", line)]
        return lines[-1] if lines else "No failure diagnostic was provided."
    return f"Process exited with code {code}; no diagnostic was provided." if type(code) is int else "No failure diagnostic was provided."


class ActivityHistory:
    """Merge an action's lifecycle into one row, scoped to its runtime attempt."""

    def __init__(self) -> None:
        self.rows: list[dict[str, Any]] = []
        self._items: dict[str, dict[str, Any]] = {}
        self._positions: dict[str, int] = {}

    def observe(self, event: Mapping[str, Any], at: str | None, *, scope: str = "0") -> dict[str, Any] | None:
        """Project known tool events, retaining details from partial updates."""
        item = event.get("item")
        if event.get("type") not in _EVENTS or not isinstance(item, Mapping):
            return
        identity = item.get("id")
        key = f"{scope}:{identity}" if isinstance(identity, str) else f"anonymous:{len(self.rows)}"
        merged = {**self._items.get(key, {}), **{k: v for k, v in item.items() if v is not None and v != ""}}
        if merged.get("type") not in _KINDS:
            return
        if merged["type"] == "agent_message":
            message = merged.get("text")
            if (event["type"] != "item.completed" or not isinstance(message, str)
                    or not message.strip() or merged.get("phase") not in {None, "commentary"}
                    or message.lstrip().startswith(("{", "[", "```")) or '"schema_version"' in message):
                return
        # Lifecycle completion wins over a stale in_progress status on partial events.
        if event["type"] in {"item.completed", "item.failed"}:
            merged["status"] = "failed" if event["type"] == "item.failed" or item.get("status") == "failed" else "completed"
        elif event["type"] == "item.started":
            merged["status"] = item.get("status") or "in_progress"
        detail = activity_detail(merged)
        row = {"kind": merged["type"], "status": merged.get("status", "in_progress"),
               "at": at, "updated_at": at, "detail": detail, "title": _title(merged, detail),
               "item_id": identity, "activity_id": key, "exit_code": merged.get("exit_code"),
               "failure_detail": _failure(merged)}
        if merged["type"] == "agent_message":
            row.update(kind="agent_update", title="Investigation update")
        if key in self._positions:
            position = self._positions[key]
            row["at"] = self.rows[position]["at"]
            self.rows[position] = row
        else:
            self._positions[key] = len(self.rows)
            self.rows.append(row)
        # Results can be large; only action fields are needed to merge later updates.
        self._items[key] = {k: v for k, v in merged.items() if k not in {"aggregated_output", "result"}}
        return row


def recover_activity(paths: list[Path], saved: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rebuild old display rows from logs without changing saved scientific evidence."""
    events: list[tuple[str, dict[str, Any]]] = []
    for path in paths:
        try:
            with path.open(encoding="utf-8", errors="replace") as stream:
                for line in stream:
                    try:
                        event = json.loads(line)
                    except ValueError:
                        continue
                    if isinstance(event, dict) and event.get("type") in _EVENTS and isinstance(event.get("item"), dict):
                        if event["item"].get("type") in _KINDS:
                            events.append((str(path), event))
        except OSError:
            continue
    if not events:
        return saved
    # Old servers kept only the last 30 events. Match timestamps backwards; do not
    # invent times for earlier actions present only in the runtime trace.
    times: dict[int, str | None] = {}
    cursor = len(events) - 1
    for row in reversed(saved):
        for index in range(cursor, -1, -1):
            event = events[index][1]
            item = event["item"]
            kind = "agent_update" if item.get("type") == "agent_message" else item.get("type")
            if kind == row.get("kind") and item.get("status", event["type"]) == row.get("status"):
                times[index] = row.get("at")
                cursor = index - 1
                break
    history = ActivityHistory()
    for index, (scope, event) in enumerate(events):
        history.observe(event, times.get(index), scope=scope)
    return history.rows
