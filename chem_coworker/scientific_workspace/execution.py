"""Recorded, bounded custom Python execution in a trusted local investigation."""

from __future__ import annotations

import json
import hashlib
import os
import subprocess
import sys
from threading import Event
from time import monotonic
from typing import Any
from uuid import uuid4

from .agent_runtime import _hidden_process_options, _stop_process
from .baseline import sha256_file, verify_baseline
from .store import canonical_bytes, InvestigationEvent, InvestigationStore


def run_python(
    store: InvestigationStore, script: str, parameters: dict[str, Any],
    evidence_refs: tuple[str, ...] = (), timeout_seconds: int = 60,
    cancel: Event | None = None,
) -> InvestigationEvent:
    """Snapshot/run Python with input/output JSON paths; record execution, not scientific truth.

    Only investigation-local scripts are accepted. This is not a security sandbox;
    the process has the investigator's OS permissions. Custom code is not replayed
    automatically. Scripts must declare external data in the investigation baseline.
    """
    root = store.root.resolve()
    path = (root / script).resolve()
    if not path.is_relative_to(root) or path.suffix != ".py":
        raise ValueError("Script must be a Python file inside this investigation")
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= 120:
        raise ValueError("timeout_seconds must be 1..120")
    if not isinstance(parameters, dict):
        raise ValueError("parameters must be a JSON object")
    if store.summary()["status"] != "active":
        raise ValueError("Investigation is stopped")
    verify_baseline(store.manifest["baseline"])
    inputs = json.loads(canonical_bytes({
        "parameters": parameters, "evidence": {ref: store.read_artifact(ref) for ref in evidence_refs},
    }))
    directory = root / "executions" / uuid4().hex
    directory.mkdir(parents=True)
    script_copy = directory / "script.py"
    script_copy.write_bytes(path.read_bytes())
    input_path, output_path = directory / "input.json", directory / "output.json"
    input_path.write_bytes(canonical_bytes(inputs))
    script_hash, input_hash = sha256_file(script_copy), sha256_file(input_path)
    payload: dict[str, Any] = {
        "schema_version": "custom_python_execution.v1", "origin": "custom_script_execution",
        "review_status": "unreviewed", "execution_directory": str(directory),
        "script": script_copy.read_text("utf-8"), "script_sha256": script_hash,
        "inputs": inputs, "input_sha256": input_hash, "interpreter": sys.executable,
        "timeout_seconds": timeout_seconds,
        "baseline_sha256": hashlib.sha256(canonical_bytes(store.manifest["baseline"])).hexdigest(),
        "limitations": ["Execution is recorded; correctness and complete input declaration are not independently verified."],
    }
    environment = os.environ.copy()
    environment["PYTHONPATH"] = store.manifest["baseline"]["repository"]
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONIOENCODING"] = "utf-8"
    options = _hidden_process_options()
    if os.name != "nt":
        options["start_new_session"] = True
    started = monotonic()
    stop = cancel or Event()
    stdout, stderr = directory / "stdout.txt", directory / "stderr.txt"
    try:
        if stop.is_set():
            raise InterruptedError("Custom execution cancelled")
        with stdout.open("wb") as out, stderr.open("wb") as err:
            process = subprocess.Popen([sys.executable, "-u", str(script_copy), str(input_path), str(output_path)],
                                       cwd=directory, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                                       env=environment, **options)
            try:
                while process.poll() is None:
                    if stop.is_set():
                        raise InterruptedError("Custom execution cancelled")
                    if monotonic() - started > timeout_seconds:
                        raise TimeoutError("Custom execution deadline exceeded")
                    if any(p.exists() and p.stat().st_size > 4_000_000 for p in (stdout, stderr, output_path)):
                        raise ValueError("Custom execution output exceeds 4 MB")
                    stop.wait(0.05)
            finally:
                _stop_process(process)
        payload["returncode"] = process.returncode
        if stop.is_set():
            raise InterruptedError("Custom execution cancelled")
        verify_baseline(store.manifest["baseline"])
        if any(p.exists() and p.stat().st_size > 4_000_000 for p in (stdout, stderr, output_path)):
            raise ValueError("Custom execution output exceeds 4 MB")
        if sha256_file(script_copy) != script_hash or sha256_file(input_path) != input_hash:
            raise ValueError("Custom script changed its recorded script or inputs")
        if process.returncode:
            raise ValueError(f"Custom script exited with code {process.returncode}")
        if output_path.stat().st_size > 4_000_000:
            raise ValueError("Custom JSON output exceeds 4 MB")
        result = json.loads(output_path.read_text("utf-8"))
        canonical_bytes(result)  # Reject non-finite JSON before recording success.
        payload["result"] = result
        payload["execution_status"] = "completed"
    except Exception as exc:
        payload["execution_status"] = "cancelled" if isinstance(exc, InterruptedError) else "timed_out" if isinstance(exc, TimeoutError) else "error"
        payload["error"] = {"type": type(exc).__name__, "message": str(exc)}
    for name, log in (("stdout", stdout), ("stderr", stderr)):
        if log.exists():
            with log.open("rb") as handle:
                payload[name] = handle.read(64000).decode("utf-8", errors="replace")
            payload[name + "_sha256"] = sha256_file(log)
            payload[name + "_truncated"] = log.stat().st_size > 64000
    payload["duration_seconds"] = round(monotonic() - started, 6)
    return store.append("custom_execution", payload, evidence_refs=evidence_refs)
