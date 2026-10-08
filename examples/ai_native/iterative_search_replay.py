"""Recorded warm-session timing control for queries already chosen by an agent.

Run via ScientificWorkspace.run_python, with parameters index_path and search_refs.
This fixed-query replay measures execution overhead, not adaptive agent quality.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from time import monotonic
from typing import Any

from condition_recommender.fragment_search import FragmentSearchSession, search_fragment_precedents


def scientific_projection(result: dict[str, Any]) -> dict[str, Any]:
    """Exclude execution telemetry added by the worker; preserve all chemistry evidence."""
    return {k: v for k, v in result.items() if k not in {"execution", "execution_status"}}


def replay(inputs: dict[str, Any]) -> dict[str, Any]:
    """Load once, rerun identical requests, and compare their complete scientific payloads."""
    parameters = inputs["parameters"]
    rows = []
    started = monotonic()
    with FragmentSearchSession(parameters["index_path"]) as session:
        session.load_library()
        load_seconds = monotonic() - started
        for ref in parameters["search_refs"]:
            payload = inputs["evidence"][ref]
            if payload.get("operation") != "search_fragment_precedents":
                raise ValueError("Replay requires saved fragment searches")
            arguments = {k: v for k, v in payload["arguments"].items()
                         if k not in {"query_variant_ref", "query_variant_id"}}
            before = monotonic()
            result = search_fragment_precedents(parameters["index_path"], **arguments, _session=session)
            rows.append({
                "search_ref": ref, "query": arguments["query"],
                "warm_seconds": round(monotonic() - before, 6),
                "cold_operation_seconds": payload["timings"]["operation_seconds"],
                "search_status": result["search_status"],
                "scientific_payload_equal": scientific_projection(result) == scientific_projection(payload["result"]),
                "original_search_status": payload["result"].get("search_status"),
                "result": result,
            })
        loads = session.library_loads
    return {"schema_version": "iterative_search_warm_replay.v1", "library_loads": loads,
            "load_seconds": round(load_seconds, 6), "total_seconds": round(monotonic() - started, 6),
            "rows": rows, "limitations": [
                "Fixed queries chosen earlier; not a live adaptive-agent comparison.",
                "Canonical matcher and ranking unchanged; private request-owned session reused experimentally.",
                "No candidate-set restriction or partial-result reuse.",
                "In-process search deadlines are cooperative; the outer recorded script has a hard timeout.",
            ]}


if __name__ == "__main__":
    request = json.loads(Path(sys.argv[1]).read_text("utf-8"))
    Path(sys.argv[2]).write_text(json.dumps(replay(request)), "utf-8")
