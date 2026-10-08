"""Deadline-limited fragment search adapter with persistent worker diagnostics."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from time import monotonic, sleep
from typing import TYPE_CHECKING, Any
from uuid import uuid4

from ..paths import REPOSITORY_ROOT

if TYPE_CHECKING:
    from .operations import ScientificOperations


def _worker_command(directory: Path) -> list[str]:
    return [sys.executable, "-m", __name__, str(directory)]


def _session_worker_command(directory: Path) -> list[str]:
    return [sys.executable, "-u", "-m", __name__, "--session", str(directory)]


class FragmentWorker:
    """One console-owned child; reuse only the index with a hard per-call deadline.

    Closing stdin or the console closes the child. Timeouts terminate the entire
    child and the next request starts fresh. There is no daemon or global cache.
    """

    def __init__(self) -> None:
        self.process: subprocess.Popen[Any] | None = None
        self.streams: list[Any] = []
        self.root: Path | None = None
        self.log_paths: dict[str, str] = {}

    def __enter__(self) -> "FragmentWorker":
        return self

    def close(self) -> None:
        """Release the owned child, including an active or unresponsive query."""
        from ..core.process_utils import stop_process_tree

        try:
            if self.process is not None:
                if self.process.poll() is None:
                    stop_process_tree(self.process)
                if self.process.stdin is not None:
                    self.process.stdin.close()
        finally:
            self.process = None
            for stream in self.streams:
                stream.close()
            self.streams.clear()

    def __exit__(self, *args: Any) -> None:
        self.close()

    def request(self, directory: Path, timeout: int) -> dict[str, Any]:
        """Submit one immutable request directory and wait for its atomic result."""
        from ..core.process_utils import hidden_process_options

        started = monotonic()
        directory = directory.resolve()
        if self.root is not None and self.root != directory.parent:
            raise ValueError("Fragment worker cannot be shared across investigations")
        if self.process is None or self.process.poll() is not None:
            self.close()
            self.root = directory.parent
            self.log_paths = {f"session_{name}": str(directory / name)
                              for name in ("stdout.log", "stderr.log")}
            self.streams = [(directory / name).open("wb") for name in ("stdout.log", "stderr.log")]
            environment = os.environ.copy()
            environment["PYTHONPATH"] = os.pathsep.join(filter(None, (
                str(REPOSITORY_ROOT), environment.get("PYTHONPATH"))))
            try:
                self.process = subprocess.Popen(
                    _session_worker_command(self.root), cwd=REPOSITORY_ROOT,
                    env=environment, stdin=subprocess.PIPE,
                    stdout=self.streams[0], stderr=self.streams[1],
                    start_new_session=os.name != "nt", **hidden_process_options(),
                )
            except BaseException:
                self.close()
                raise
        try:
            assert self.process.stdin is not None
            self.process.stdin.write((json.dumps(str(directory)) + "\n").encode("utf-8"))
            self.process.stdin.flush()
            while not (directory / "result.json").exists():
                if self.process.poll() is not None:
                    raise RuntimeError("Fragment session worker exited before returning evidence")
                if monotonic() - started >= timeout:
                    raise subprocess.TimeoutExpired(self.process.args, timeout)
                sleep(0.02)
            result = json.loads((directory / "result.json").read_text("utf-8"))
            if result.get("execution_status") not in {"completed", "error", "timed_out"}:
                raise RuntimeError("Invalid fragment session result")
            if result["execution_status"] == "timed_out":
                self.close()
            return result
        except BaseException:
            self.close()
            raise


def run_fragment_search(operations: ScientificOperations, arguments: dict[str, Any], *, automatic: bool = False) -> dict[str, Any]:
    """Run one prepared-index query in a killable child; never build an index."""
    from ..core.process_utils import hidden_process_options, stop_process_tree

    timeout = arguments["timeout_seconds"]
    maximum = 120 if automatic else 30
    if type(timeout) is not int or not 1 <= timeout <= maximum:
        raise ValueError(f"timeout_seconds must be an integer between 1 and {maximum}")
    path = operations._path("fragment_index")
    directory = operations.store.root / "diagnostics" / "fragment_search" / uuid4().hex
    directory.mkdir(parents=True)
    request = {"index_path": str(path), **arguments}
    if automatic:
        request["operation"] = "find_synthesis_precedents"
    (directory / "request.json").write_text(json.dumps(request), "utf-8")
    repository = str(REPOSITORY_ROOT)
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, (repository, environment.get("PYTHONPATH"))))
    started = monotonic()
    result: dict[str, Any] = {}
    process = None
    try:
        worker = None if automatic else operations.fragment_worker
        if worker is not None:
            result = worker.request(directory, timeout)
        else:
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
    if not automatic and operations.fragment_worker is not None:
        result["execution"]["diagnostics"].update(operations.fragment_worker.log_paths)
    return result


def _worker(directory: Path, session: Any = None, session_id: str | None = None) -> None:
    from condition_recommender.fragment_search import search_fragment_precedents

    def progress(event: dict[str, Any]) -> None:
        with (directory / "stages.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event) + "\n")

    try:
        request = json.loads((directory / "request.json").read_text("utf-8"))
        operation = request.pop("operation", "search_fragment_precedents")
        if operation == "find_synthesis_precedents":
            from condition_recommender.precedent_discovery import find_synthesis_precedents
            result = find_synthesis_precedents(**request, progress=progress)
        elif operation == "search_fragment_precedents":
            reused = session is not None and session.library is not None
            result = search_fragment_precedents(**request, progress=progress, _session=session)
            if session is not None:
                result["execution"].update(session_id=session_id, library_reused=reused,
                                           session_library_loads=session.library_loads)
        else:
            raise ValueError("Unsupported fragment worker operation")
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


def _session_worker(root: Path) -> None:
    from condition_recommender.fragment_search import FragmentSearchSession

    session = None
    session_id = uuid4().hex
    try:
        for line in sys.stdin:
            directory = Path(json.loads(line)).resolve()
            if directory.parent != root.resolve():
                raise ValueError("Request is outside the owning investigation")
            request = json.loads((directory / "request.json").read_text("utf-8"))
            if session is None:
                session = FragmentSearchSession(request["index_path"])
            _worker(directory, session, session_id)
            # Agent edits need full searches, never an implicit narrowed candidate set.
            session.candidates.clear()
    finally:
        if session is not None:
            session.__exit__(None, None, None)


if __name__ == "__main__":
    if sys.argv[1] == "--session":
        _session_worker(Path(sys.argv[2]))
    else:
        _worker(Path(sys.argv[1]))
