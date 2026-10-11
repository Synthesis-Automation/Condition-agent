"""Isolated JSON-file worker, executable under RetroChimera's Python 3.9.

No project imports: this script runs in the separately provisioned ML environment.
The parent owns the process tree and deadline, including loading and hashing.
"""

from __future__ import annotations

import hashlib
from importlib.metadata import version
import json
import multiprocessing
import os
from pathlib import Path
import sys
from time import monotonic
import traceback
from typing import Any


def fingerprint(path: Path) -> str:
    """Hash a dependency without reading a checkpoint into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def token(path: Path) -> tuple[int, int]:
    """Detect dependency changes during a worker session."""
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def main() -> None:
    """Load once and serve serial requests until the parent closes stdin."""
    config = json.loads(Path(sys.argv[1]).read_text("utf-8"))
    model = None
    snapshots = {}
    identity = {}
    load_seconds = 0.0
    for line in sys.stdin:
        directory = Path(json.loads(line))
        started = monotonic()
        reused = model is not None
        try:
            request = json.loads((directory / "request.json").read_text("utf-8"))
            if model is None:
                for entry in config["pinned_files"]:
                    path = Path(entry["path"])
                    before = token(path)
                    if fingerprint(path) != entry["sha256"] or token(path) != before:
                        raise ValueError(f"Pinned supplier dependency changed: {path}")
                    snapshots[path] = before
                os.environ["OMP_NUM_THREADS"] = str(config["num_threads"])
                os.environ["MKL_NUM_THREADS"] = str(config["num_threads"])
                sys.path.insert(0, config["checkout"])
                import torch
                from syntheseus import Molecule
                from retrochimera import RetroChimeraModel, inference
                from retrochimera.chem.rules import RuleBasedRetrosynthesizer

                torch.set_num_threads(config["num_threads"])
                torch.manual_seed(0)
                if config["device"] == "cuda" and not torch.cuda.is_available():
                    raise RuntimeError("Configured CUDA device is unavailable")
                model_dir = Path(config["model_dir"])
                components = json.loads((model_dir / "models.json").read_text("utf-8"))
                kwargs = {
                    name: {"num_processes": config["num_processes"]}
                    for name, (class_name, _) in components.items()
                    if issubclass(getattr(inference, class_name), RuleBasedRetrosynthesizer)
                }
                model = RetroChimeraModel(
                    model_dir=model_dir, device=config["device"],
                    call_submodels_in_parallel=False, use_cache=False, **kwargs,
                )
                identity = {
                    "python": sys.version, "versions": {name: version(name) for name in (
                        "torch", "rdkit", "syntheseus", "retrochimera")},
                    "components": {name: spec[0] for name, spec in components.items()},
                    "config_sha256": fingerprint(Path(sys.argv[1])),
                    "device": config["device"], "seed": 0,
                }
                load_seconds = monotonic() - started
            if any(token(path) != previous for path, previous in snapshots.items()):
                raise ValueError("Supplier files changed; prepare a new configuration and investigation")
            model.consensus_only = request["consensus_only"]
            inference_started = monotonic()
            predictions = model([Molecule(request["target_smiles"])], num_results=request["num_results"])[0]
            rows = [{
                "product": p.product.smiles, "reactants": sorted(m.smiles for m in p.reactants),
                "score": p.metadata.get("score"), "probability": p.metadata.get("probability"),
                "individual_ranks": p.metadata.get("individual_ranks", {}),
            } for p in predictions]
            result = {"execution_status": "completed", "predictions": rows, "identity": identity,
                      "execution": {"model_reused": reused, "model_load_seconds": 0 if reused else load_seconds,
                                    "inference_seconds": monotonic() - inference_started}}
        except Exception as exc:
            traceback.print_exc()
            result = {"execution_status": "error", "predictions": [],
                      "error": {"type": type(exc).__name__, "message": str(exc)}}
        result.setdefault("execution", {}).update(worker_pid=os.getpid(), elapsed_seconds=monotonic() - started)
        temporary = directory / "result.tmp"
        temporary.write_text(json.dumps(result, allow_nan=False), "utf-8")
        temporary.replace(directory / "result.json")
        if result["execution_status"] != "completed":
            break


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
