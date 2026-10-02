"""Identify the executable scientific environment without changing its data."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import sqlite3
import subprocess
from contextlib import closing
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

from ..paths import REPOSITORY_ROOT

PACKAGE_ROOTS = (
    "reactive_taxonomy", "condition_registry", "condition_recommender",
    "core_retrosynthesis", "forward_synthesis", "cas_tools", "chem_coworker",
)

BASELINE_SCHEMA = "scientific_baseline.v2"
LEGACY_BASELINE_SCHEMA = "scientific_baseline.v1"
WORKSPACE_PREFIX = "chem_coworker/scientific_workspace/"

# Only these application responsibilities may change without replacing the
# scientific baseline. Unknown modules, including new tool adapters, are
# scientific by default. Evidence validators and execution stay scientific.
APPLICATION_MODULE_LAYERS: Mapping[str, str] = MappingProxyType({
    "runtime/__init__.py": "runtime",
    "agent_context/__init__.py": "guidance",
    "runtime/activity.py": "runtime",
    "runtime/agent_runtime.py": "runtime",
    "answers/answer_handoff.py": "runtime",
    "runtime/conversation.py": "runtime",
    "runtime/research_profiles.py": "runtime",
    "runtime/runtime_environment.py": "runtime",
    "agent_context/context.py": "guidance",
    "agent_context/learning.py": "guidance",
    "agent_context/prompts.py": "guidance",
})
APPLICATION_DIRECTORY_LAYERS: Mapping[str, str] = MappingProxyType({
    "task_playbooks": "guidance", "agent_instructions": "guidance", "presentation": "presentation",
})


def application_layer(path: str) -> str | None:
    """Classify explicitly owned application files; everything else is scientific."""
    if not path.startswith(WORKSPACE_PREFIX):
        return None
    relative = path[len(WORKSPACE_PREFIX):]
    if relative in APPLICATION_MODULE_LAYERS:
        return APPLICATION_MODULE_LAYERS[relative]
    if Path(relative).suffix != ".md":
        return None
    return APPLICATION_DIRECTORY_LAYERS.get(relative.split("/", 1)[0])


def sha256_file(path: Path) -> str:
    """Hash an artifact in bounded memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def legacy_code_manifest(repository: Path) -> dict[str, str]:
    """Keep v1 strict for both original and renamed task resource directories."""
    workspace_directory = repository / "chem_coworker" / "scientific_workspace"
    task_directories = {workspace_directory / "guides", workspace_directory / "task_playbooks"}
    files = {}
    for name in PACKAGE_ROOTS:
        for path in sorted((repository / name).rglob("*")):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.suffix == ".py" or (
                "definitions" in path.parts and path.suffix in {".json", ".jsonl"}
            ) or (
                path.parent in task_directories
                and path.suffix == ".md"
            ):
                files[path.relative_to(repository).as_posix()] = sha256_file(path)
    return files


def code_manifest(repository: Path) -> dict[str, str]:
    """Hash scientific code and definitions, conservatively including new adapters.

    Runtime configuration, task guidance, and presentation are separately recorded
    application context. Their explicitly owned files do not govern scientific
    replay; all other Python modules and workspace JSON contracts do.
    """
    files = {}
    for name in PACKAGE_ROOTS:
        for path in sorted((repository / name).rglob("*")):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            relative = path.relative_to(repository).as_posix()
            if application_layer(relative) is not None:
                continue
            if path.suffix == ".py" or (
                path.suffix in {".json", ".jsonl"}
                and ("definitions" in path.parts or relative.startswith(WORKSPACE_PREFIX))
            ):
                files[relative] = sha256_file(path)
    return files


def scientific_identity(baseline: Mapping[str, Any]) -> str:
    """Identify scientific inputs independently of guidance and display settings."""
    from .store import canonical_bytes

    identity = {
        "schema_version": baseline.get("schema_version", BASELINE_SCHEMA),
        "code_files": baseline["code_files"],
        "environment": baseline["environment"],
        "taxonomy_versions": baseline.get("taxonomy_versions", {}),
        "registry_versions": baseline.get("registry_versions", {}),
        "artifacts": {
            name: {key: value for key, value in artifact.items() if key != "mtime_ns"}
            for name, artifact in baseline["artifacts"].items()
        },
    }
    return "sha256:" + hashlib.sha256(canonical_bytes(identity)).hexdigest()


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
    for package in ("rdkit", "numpy", "rdchiral", "pypdf"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed"
    return versions


def capture_baseline(
    repository: Path, artifacts: Mapping[str, Path],
) -> dict[str, Any]:
    """Record a development snapshot, registry audit, and selected data inputs."""
    from condition_registry import (
        condition_registry_definition_versions,
        validate_registry,
    )
    from reactive_taxonomy import reaction_signature_definition_versions

    repository = repository.resolve()
    if repository != REPOSITORY_ROOT:
        raise ValueError("Repository must match the imported scientific packages")
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repository, text=True,
        capture_output=True, check=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"], cwd=repository, text=True,
        capture_output=True, check=True,
    ).stdout.splitlines()
    paths = dict(artifacts)
    if "weak_label_records" in paths:
        from condition_recommender import weak_label_recipe_catalog_path

        catalog = weak_label_recipe_catalog_path(paths["weak_label_records"]).resolve()
        configured = paths.get("weak_label_recipe_catalog", catalog).resolve()
        if configured != catalog:
            raise ValueError("weak_label_recipe_catalog must be the catalog beside weak_label_records")
        paths["weak_label_recipe_catalog"] = catalog
    identities = {name: artifact_identity(path) for name, path in paths.items()}
    for value in identities.values():
        path = Path(value["path"])
        if value["status"] == "present" and path.suffix == ".sqlite":
            with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as connection:
                tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                if "metadata" in tables:
                    columns = {row[1] for row in connection.execute("PRAGMA table_info(metadata)")}
                    if "payload" in columns:
                        value["metadata"] = [json.loads(row[0]) for row in connection.execute("SELECT payload FROM metadata")]
    baseline = {
        "schema_version": BASELINE_SCHEMA,
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
    from ..agent_context.learning import build_learning_context

    baseline["learning_context"] = build_learning_context(baseline)
    baseline["scientific_identity"] = scientific_identity(baseline)
    return baseline


def verify_baseline(baseline: Mapping[str, Any], *, full_hash: bool = False) -> None:
    """Reject code, dependency, or input changes before continuing an investigation.

    Normal calls inspect data stat identity; replay also rehashes large inputs.
    No changed artifact is automatically accepted as a new scientific baseline.
    """
    schema = baseline.get("schema_version", BASELINE_SCHEMA)
    if schema not in {BASELINE_SCHEMA, LEGACY_BASELINE_SCHEMA}:
        raise ValueError("Unsupported scientific baseline schema")
    if ("scientific_identity" in baseline
            and scientific_identity(baseline) != baseline["scientific_identity"]):
        raise ValueError("Scientific baseline identity mismatch")
    manifest = legacy_code_manifest if schema == LEGACY_BASELINE_SCHEMA else code_manifest
    if manifest(Path(baseline["repository"])) != baseline["code_files"]:
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
