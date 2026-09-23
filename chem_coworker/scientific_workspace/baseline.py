"""Identify the executable scientific environment without changing its data."""

from __future__ import annotations

import hashlib
import importlib.metadata
import platform
import json
import sqlite3
import subprocess
from contextlib import closing
from pathlib import Path
from typing import Any, Mapping


PACKAGE_ROOTS = (
    "reactive_taxonomy", "condition_registry", "condition_recommender",
    "core_retrosynthesis", "cas_tools", "chem_coworker",
)


def sha256_file(path: Path) -> str:
    """Hash an artifact in bounded memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def code_manifest(repository: Path) -> dict[str, str]:
    """Hash code and definitions, including uncommitted implementation changes."""
    files = {}
    for name in PACKAGE_ROOTS:
        for path in sorted((repository / name).rglob("*")):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.suffix == ".py" or (
                "definitions" in path.parts and path.suffix in {".json", ".jsonl"}
            ):
                files[path.relative_to(repository).as_posix()] = sha256_file(path)
    return files


def artifact_identity(path: Path) -> dict[str, Any]:
    """Capture full content identity; report missing inputs explicitly."""
    path = path.expanduser().resolve()
    if not path.is_file():
        return {"path": str(path), "status": "missing"}
    before = path.stat()
    digest = sha256_file(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Artifact changed during fingerprinting: {path}")
    return {
        "path": str(path), "status": "present", "sha256": digest,
        "size_bytes": after.st_size, "mtime_ns": after.st_mtime_ns,
    }


def environment_versions() -> dict[str, str]:
    """Record runtime versions relevant to deterministic calculations."""
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for package in ("rdkit", "numpy", "rdchiral"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed"
    return versions


def capture_baseline(
    repository: Path, artifacts: Mapping[str, Path],
) -> dict[str, Any]:
    """Record a development snapshot, registry audit, and selected data inputs."""
    from condition_registry import condition_registry_definition_versions, validate_registry
    from reactive_taxonomy import reaction_signature_definition_versions

    repository = repository.resolve()
    if repository != Path(__file__).resolve().parents[2]:
        raise ValueError("Repository must match the imported scientific packages")
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repository, text=True,
        capture_output=True, check=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"], cwd=repository, text=True,
        capture_output=True, check=True,
    ).stdout.splitlines()
    identities = {name: artifact_identity(path) for name, path in artifacts.items()}
    for value in identities.values():
        path = Path(value["path"])
        if value["status"] == "present" and path.suffix == ".sqlite":
            with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as connection:
                tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                if "metadata" in tables:
                    columns = {row[1] for row in connection.execute("PRAGMA table_info(metadata)")}
                    if "payload" in columns:
                        value["metadata"] = [json.loads(row[0]) for row in connection.execute("SELECT payload FROM metadata")]
    return {
        "schema_version": "scientific_baseline.v1",
        "repository": str(repository), "git_revision": revision,
        "git_status": dirty, "code_files": code_manifest(repository),
        "environment": environment_versions(),
        "taxonomy_versions": reaction_signature_definition_versions(),
        "registry_versions": condition_registry_definition_versions(),
        "registry_validation": validate_registry(),
        "artifacts": identities,
        "validation_status": "development_snapshot_not_release_validated",
        "evaluation_partition": "development_only_not_an_untouched_evaluation",
    }


def verify_baseline(baseline: Mapping[str, Any], *, full_hash: bool = False) -> None:
    """Reject code, dependency, or input changes before continuing an investigation.

    Normal calls inspect data stat identity; replay also rehashes large inputs.
    No changed artifact is automatically accepted as a new scientific baseline.
    """
    if code_manifest(Path(baseline["repository"])) != baseline["code_files"]:
        raise ValueError("Scientific code or definitions changed; start a new investigation")
    if environment_versions() != baseline["environment"]:
        raise ValueError("Scientific runtime changed; start a new investigation")
    for value in baseline["artifacts"].values():
        path = Path(value["path"])
        if value["status"] == "missing":
            if path.exists():
                raise ValueError(f"Previously missing artifact appeared: {path}")
            continue
        if not path.is_file():
            raise ValueError(f"Baseline artifact is missing: {path}")
        stat = path.stat()
        if (stat.st_size, stat.st_mtime_ns) != (value["size_bytes"], value["mtime_ns"]):
            raise ValueError(f"Baseline artifact changed: {path}")
        if full_hash and sha256_file(path) != value["sha256"]:
            raise ValueError(f"Baseline artifact content changed: {path}")
