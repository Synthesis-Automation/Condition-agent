"""Optional, session-owned RetroChimera process and recorded workspace adapter."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
from threading import Lock
from time import monotonic, sleep
from typing import Any, TYPE_CHECKING
from uuid import uuid4

from core_retrosynthesis.supplier_predictions import canonical_supplier_smiles, normalize_supplier_predictions
from ..core.baseline import sha256_file
from ..core.process_utils import hidden_process_options, stop_process_tree
from .retrochimera_config import load_config

if TYPE_CHECKING:
    from .operations import ScientificOperations


class RetroChimeraWorker:
    """Own a serial persistent worker; deadline/error closes its entire process tree."""

    def __init__(self) -> None:
        self.process: subprocess.Popen[Any] | None = None
        self.streams: list[Any] = []
        self.session_key: tuple[str, str] | None = None
        self.log_paths: list[str] = []
        self.lock = Lock()

    def __enter__(self) -> RetroChimeraWorker:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def close(self) -> None:
        """Stop only this worker and its rule-application children."""
        try:
            if self.process is not None:
                stop_process_tree(self.process)
                if self.process.stdin:
                    self.process.stdin.close()
        finally:
            self.process = None
            for stream in self.streams:
                stream.close()
            self.streams.clear()
            self.session_key = None

    def request(self, config_path: Path, directory: Path, timeout: int) -> dict[str, Any]:
        """Submit a request with a deadline covering startup, hashing and inference."""
        if not self.lock.acquire(blocking=False):
            raise RuntimeError("RetroChimera worker already has an active request")
        started = monotonic()
        try:
            config = load_config(config_path)
            if config.get("device") not in {"cpu", "cuda"}:
                raise ValueError("Configure device cpu or cuda")
            for key in ("num_threads", "num_processes"):
                if type(config.get(key)) is not int or not 1 <= config[key] <= 32:
                    raise ValueError(f"{key} must be an integer in 1..32")
            if not config.get("pinned_files"):
                raise ValueError("Supplier configuration must pin its model and source files")
            key = (str(directory.parent.resolve()), sha256_file(config_path))
            if self.session_key != key or self.process is None or self.process.poll() is not None:
                self.close()
                self.session_key = key
                self.log_paths = [str(directory / f"worker_{name}.log") for name in ("stdout", "stderr")]
                self.streams = [Path(path).open("wb") for path in self.log_paths]
                environment = os.environ.copy()
                environment["PYTHONPATH"] = config["checkout"]
                environment["PYTHONNOUSERSITE"] = "1"
                self.process = subprocess.Popen(
                    [config["python_executable"], "-u", str(Path(__file__).with_name("retrochimera_worker.py")), str(config_path)],
                    cwd=config["checkout"], env=environment, stdin=subprocess.PIPE,
                    stdout=self.streams[0], stderr=self.streams[1],
                    start_new_session=os.name != "nt", **hidden_process_options(),
                )
            assert self.process.stdin is not None
            self.process.stdin.write((json.dumps(str(directory.resolve())) + "\n").encode("utf-8"))
            self.process.stdin.flush()
            while not (directory / "result.json").is_file():
                if self.process.poll() is not None:
                    raise RuntimeError("RetroChimera worker exited without a result; inspect worker logs")
                if monotonic() - started >= timeout:
                    raise subprocess.TimeoutExpired(self.process.args, timeout)
                sleep(0.02)
            result = json.loads((directory / "result.json").read_text("utf-8"))
            if result.get("execution_status") not in {"completed", "error"}:
                raise ValueError("Invalid RetroChimera worker response")
            if result["execution_status"] != "completed":
                self.close()
        except subprocess.TimeoutExpired:
            self.close()
            result = {"execution_status": "timed_out", "predictions": [],
                      "error": {"type": "TimeoutExpired", "message": f"RetroChimera exceeded {timeout}s including startup"}}
        except KeyboardInterrupt:
            self.close()
            raise
        except Exception as exc:
            self.close()
            result = {"execution_status": "error", "predictions": [],
                      "error": {"type": type(exc).__name__, "message": str(exc)}}
        finally:
            self.lock.release()
        result.setdefault("execution", {}).update(
            elapsed_seconds=monotonic() - started, worker_logs=list(self.log_paths),
        )
        return result


def predict_retrochimera(
    operations: ScientificOperations, target_smiles: str, num_results: int,
    timeout_seconds: int, consensus_only: bool,
) -> dict[str, Any]:
    """Record model proposals separately from the existing operator supplier."""
    if type(num_results) is not int or not 1 <= num_results <= 50:
        raise ValueError("num_results must be an integer in 1..50; prefer 5..10")
    if type(timeout_seconds) is not int or not 1 <= timeout_seconds <= 600:
        raise ValueError("timeout_seconds must be an integer in 1..600")
    if type(consensus_only) is not bool:
        raise ValueError("consensus_only must be boolean")
    target = canonical_supplier_smiles(target_smiles)
    config_path = operations._path("retrochimera_config")
    directory = operations.store.root / "diagnostics" / "retrochimera" / uuid4().hex
    directory.mkdir(parents=True)
    request = {"target_smiles": target, "num_results": num_results, "consensus_only": consensus_only}
    (directory / "request.json").write_text(json.dumps(request), "utf-8")
    worker = operations.retrochimera_worker
    if worker is None:
        with RetroChimeraWorker() as temporary:
            result = temporary.request(config_path, directory, timeout_seconds)
    else:
        result = worker.request(config_path, directory, timeout_seconds)
    normalized = normalize_supplier_predictions(target, result.get("predictions", []))
    return {
        **normalized, "execution_status": result["execution_status"], "request": request,
        "identity": result.get("identity", {}), "execution": result.get("execution", {}),
        "config_sha256": sha256_file(config_path),
        **({"error": result["error"]} if "error" in result else {}),
    }
