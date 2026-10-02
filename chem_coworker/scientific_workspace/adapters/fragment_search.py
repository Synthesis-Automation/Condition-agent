"""Deadline-limited fragment search adapter with persistent worker diagnostics."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from time import monotonic
from typing import TYPE_CHECKING, Any
from uuid import uuid4

from ..paths import REPOSITORY_ROOT

if TYPE_CHECKING:
    from .operations import ScientificOperations


def _worker_command(directory: Path) -> list[str]:
    return [sys.executable, "-m", __name__, str(directory)]


def run_fragment_search(operations: ScientificOperations, arguments: dict[str, Any]) -> dict[str, Any]:
    """Run one prepared-index query in a killable child; never build an index."""
    from ..core.process_utils import hidden_process_options, stop_process_tree

    timeout = arguments["timeout_seconds"]
    if type(timeout) is not int or not 1 <= timeout <= 30:
        raise ValueError("timeout_seconds must be an integer between 1 and 30")
    path = operations._path("fragment_index")
    directory = operations.store.root / "diagnostics" / "fragment_search" / uuid4().hex
    directory.mkdir(parents=True)
    request = {"index_path": str(path), **arguments}
    (directory / "request.json").write_text(json.dumps(request), "utf-8")
    repository = str(REPOSITORY_ROOT)
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, (repository, environment.get("PYTHONPATH"))))
    started = monotonic()
    result: dict[str, Any] = {}
    process = None
    try:
        with (directory / "stdout.log").open("wb") as stdout, (directory / "stderr.log").open("wb") as stderr:
            process = subprocess.Popen(_worker_command(directory),
                                       cwd=repository, env=environment, stdout=stdout, stderr=stderr,
                                       start_new_session=os.name != "nt", **hidden_process_options())
            process.wait(timeout=max(0, timeout - (monotonic() - started)))
        result = json.loads((directory / "result.json").read_text("utf-8"))
        if process.returncode != 0 or result.get("execution_status") not in {"completed", "error", "timed_out"}:
            raise RuntimeError("Fragment worker returned an invalid result")
    except subprocess.TimeoutExpired:
        result = {"execution_status": "timed_out", "search_status": "partial",
                  "stop_reason": "worker_deadline", "hits": [],
                  "counts": {"products": {"value": None, "precision": "unknown"}},
                  "error": {"type": "TimeoutExpired", "message": f"Fragment search exceeded {timeout} seconds; no absence claim is supported."}}
    except KeyboardInterrupt:
        result = {"execution_status": "cancelled", "error": {"type": "KeyboardInterrupt", "message": "Fragment search interrupted"}}
    finally:
        if process is not None and process.poll() is None:
            stop_process_tree(process)
    result.setdefault("execution", {}).update({
        "worker_elapsed_seconds": round(monotonic() - started, 6),
        "diagnostics": {p.name: str(p) for p in sorted(directory.iterdir()) if p.is_file()},
    })
    return result


def _worker(directory: Path) -> None:
    from condition_recommender.fragment_search import search_fragment_precedents

    def progress(event: dict[str, Any]) -> None:
        with (directory / "stages.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event) + "\n")

    try:
        request = json.loads((directory / "request.json").read_text("utf-8"))
        result = search_fragment_precedents(**request, progress=progress)
        timed_out = "deadline" in str(result.get("stop_reason"))
        result["execution_status"] = "timed_out" if timed_out else "completed"
        if timed_out:
            result["error"] = {"type": "TimeoutExpired", "message": "Fragment search returned partial evidence at its deadline"}
    except Exception as exc:
        result = {"execution_status": "error", "error": {"type": type(exc).__name__, "message": str(exc)}}
        progress({"stage": "error", **result["error"]})
    temporary = directory / "result.tmp"
    temporary.write_text(json.dumps(result), "utf-8")
    temporary.replace(directory / "result.json")


if __name__ == "__main__":
    _worker(Path(sys.argv[1]))
