"""Line-oriented agent console with one owned, deadline-limited fragment worker.

Run ``python -u -m chem_coworker.scientific_workspace.console WORKSPACE`` in a
persistent terminal session. Send one JSON request per line, inspect its response,
then choose the next request. EOF, quit, interruption or the call cap closes the
child. This adapter schedules no chemistry and selects no fragments.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from .adapters.fragment_search import FragmentWorker
from .adapters.operations import ScientificOperations
from .core.store import InvestigationStore
from .workspace import ScientificWorkspace


class ScientificConsole:
    """Translate single agent requests into canonical recorded workspace calls."""

    def __init__(self, workspace: ScientificWorkspace) -> None:
        self.workspace = workspace

    def dispatch(self, request: dict[str, Any]) -> dict[str, Any]:
        """Run, inspect, or describe; preserve the reason and parent evidence."""
        w = self.workspace
        action = request.get("action", "run")
        if action == "help":
            return {"help": w.help(request["names"])}
        refs = tuple(request.get("evidence_refs", []))
        if request.get("parent_ref"):
            refs += (request["parent_ref"],)
        for ref in refs:
            w.store.read_artifact(ref)
        if action == "inspect":
            ref = request["source_ref"]
            view = w.inspect_artifact(ref, tuple(request.get("path", ("result",))),
                                      offset=request.get("offset", 0), limit=request.get("limit", 3))
            event = w.store.append("console_inspection", {
                "schema_version": "scientific_console_inspection.v1",
                "source_ref": ref, "path": request.get("path", ["result"]),
                "offset": request.get("offset", 0), "limit": request.get("limit", 3),
                "reason": request.get("reason"), "review_status": "opened_not_adjudicated",
            }, evidence_refs=(ref, *refs))
            return {"inspection_ref": event.artifact_ref, "view": view}
        if action == "search_prepared":
            ref = request["query_ref"]
            payload = w.store.read_artifact(ref)
            if (payload.get("operation") != "propose_fragment_queries"
                    or payload.get("execution_status") != "completed"):
                raise ValueError("query_ref must identify a completed query preview")
            preview = payload["result"]
            variant_id = request.get("variant_id")
            if variant_id:
                variants = [v for v in preview["variants"] if v["variant_id"] == variant_id]
                if len(variants) != 1:
                    raise ValueError("variant_id is absent from the saved preview")
                selected = variants[0]
            else:
                selected = {**preview["parent_query"], "query": preview["parent_query"]["expression"]}
            arguments = {k: selected[k] for k in ("query", "query_format", "topology")}
            arguments.update(target_smiles=preview["target_smiles"], limit=request.get("limit", 5),
                             timeout_seconds=request.get("timeout_seconds", 30))
            if variant_id:
                arguments.update(query_variant_ref=ref, query_variant_id=variant_id)
            operation = "search_fragment_precedents"
            refs += (ref,)
        elif action == "run":
            operation, arguments = request["operation"], request.get("arguments", {})
        else:
            raise ValueError("Use run, search_prepared, inspect, help, or quit")
        reason = str(request.get("reason", "")).strip()
        if operation in {"search_fragment_precedents", "propose_fragment_queries"} and not reason:
            raise ValueError("State the chemical search question or query-revision reason")
        if reason:
            w.store.note("decision", reason, evidence_refs=refs)
        event = w.run(operation, arguments, evidence_refs=refs)
        response = w.call_summary(event)
        payload = w.store.read_artifact(event.artifact_ref)
        execution = (payload.get("result") or {}).get("execution", {})
        response["search_execution"] = {k: execution[k] for k in (
            "session_id", "library_reused", "session_library_loads", "worker_elapsed_seconds") if k in execution}
        response["timings"] = payload.get("timings", {})
        return response


def main() -> None:
    """Serve bounded JSONL requests over stdin without a listening socket or daemon."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path)
    parser.add_argument("--max-requests", type=int, default=40)
    args = parser.parse_args()
    if not 1 <= args.max_requests <= 200:
        parser.error("max-requests must be in 1..200")
    store = InvestigationStore(args.workspace)
    with FragmentWorker() as worker:
        operations = ScientificOperations(store, fragment_worker=worker)
        console = ScientificConsole(ScientificWorkspace(args.workspace, operations=operations))
        print(json.dumps({"status": "ready", "protocol": "scientific_console.v1",
                          "max_requests": args.max_requests}), flush=True)
        for _ in range(args.max_requests):
            line = sys.stdin.readline(65537)
            if not line:
                break
            if len(line) > 65536:
                print(json.dumps({"error": "Request exceeds 64 KiB; console closed"}), flush=True)
                break
            try:
                request = json.loads(line)
                if not isinstance(request, dict):
                    raise ValueError("Each request must be a JSON object")
                if request.get("action") == "quit":
                    break
                result = console.dispatch(request)
                print(json.dumps(result), flush=True)
            except Exception as exc:
                print(json.dumps({"error": {"type": type(exc).__name__, "message": str(exc)}}), flush=True)
        print(json.dumps({"status": "closed"}), flush=True)


if __name__ == "__main__":
    main()
