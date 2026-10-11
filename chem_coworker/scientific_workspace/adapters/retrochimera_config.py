"""Prepare a pinned local RetroChimera installation; never install or download."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

from ..core.baseline import sha256_file


def load_config(path: Path) -> dict[str, Any]:
    """Validate that every model asset and supplier module has a pinned identity."""
    config = json.loads(path.read_text("utf-8"))
    if config.get("schema_version") != "retrochimera_config.v1":
        raise ValueError("Unsupported RetroChimera configuration schema")
    checkout = Path(config["checkout"])
    model_dir = Path(config["model_dir"])
    python = Path(config["python_executable"])
    if not all(p.is_absolute() for p in (checkout, model_dir, python)):
        raise ValueError("RetroChimera configuration paths must be absolute")
    if not (checkout / "retrochimera" / "__init__.py").is_file():
        raise FileNotFoundError("RetroChimera checkout is unavailable")
    if not (model_dir / "models.json").is_file() or not python.is_file():
        raise FileNotFoundError("RetroChimera model or interpreter is unavailable")
    required = {python.resolve(), *(checkout / "retrochimera").rglob("*.py"),
                *(p for p in model_dir.rglob("*") if p.is_file())}
    entries = config.get("pinned_files", [])
    paths = []
    for entry in entries:
        dependency = Path(entry["path"])
        if not dependency.is_absolute() or not re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]):
            raise ValueError("Invalid pinned supplier dependency")
        paths.append(dependency.resolve())
    if len(set(paths)) != len(paths) or set(paths) != required:
        raise ValueError("Pinned supplier files do not match the complete installation/model assets")
    return config


def prepare_config(checkout: Path, python: Path, model_dir: Path) -> dict[str, Any]:
    """Fingerprint model assets and supplier source for a trusted local worker."""
    checkout, python, model_dir = (p.resolve() for p in (checkout, python, model_dir))
    if not python.is_file() or not (model_dir / "models.json").is_file():
        raise FileNotFoundError("Supply an existing Python executable and extracted checkpoint")
    source = checkout / "retrochimera"
    if not (source / "__init__.py").is_file():
        raise FileNotFoundError("checkout must contain the RetroChimera Python package")
    files = sorted(set([python, *source.rglob("*.py"), *[p for p in model_dir.rglob("*") if p.is_file()]]))
    return {
        "schema_version": "retrochimera_config.v1", "checkout": str(checkout),
        "python_executable": str(python), "model_dir": str(model_dir),
        "device": "cpu", "num_processes": 2, "num_threads": 4,
        "pinned_files": [{"path": str(p), "sha256": sha256_file(p)} for p in files],
    }


def main() -> None:
    """Write an explicit configuration to be registered as a workspace artifact."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    config = prepare_config(args.checkout, args.python, args.model_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2)
    print(args.output)


if __name__ == "__main__":
    main()
