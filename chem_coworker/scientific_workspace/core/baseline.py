"""Identify the executable scientific environment without changing its data."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import sqlite3
import subprocess
import sys
import sysconfig
from contextlib import closing
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Mapping

from ..paths import REPOSITORY_ROOT
from .fingerprint_cache import FingerprintCache, stat_token

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
    "runtime/workspace_modes.py": "runtime",
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


def sha256_file(path: Path, *, on_progress: Callable[[int], None] | None = None) -> str:
    """Hash in bounded memory; the optional byte callback may abort the read."""
    digest = hashlib.sha256()
    bytes_read = 0
    if on_progress is not None:
        on_progress(bytes_read)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            bytes_read += len(chunk)
            if on_progress is not None:
                on_progress(bytes_read)
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


def artifact_identity(
    path: Path, *, on_progress: Callable[[int], None] | None = None,
    cache: FingerprintCache | None = None,
    on_reuse: Callable[[], None] | None = None,
    expected_sha256: str | None = None,
) -> dict[str, Any]:
    """Capture content identity, optionally reusing a stat-validated local hash."""
    path = path.expanduser().resolve()
    if not path.is_file():
        return {"path": str(path), "status": "missing"}
    before = path.stat()
    digest = cache.get(path, before) if cache is not None else None
    if expected_sha256 is not None and digest != expected_sha256:
        digest = None
    reused = digest is not None
    if digest is None:
        digest = sha256_file(path, on_progress=on_progress)
    elif on_reuse is not None:
        on_reuse()
    after = path.stat()
    if stat_token(before) != stat_token(after):
        raise ValueError(f"Artifact changed during fingerprinting: {path}")
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError(f"Processed artifact checksum mismatch: {path}")
    if cache is not None and not reused:
        cache.put(path, after, digest)
    return {
        "path": str(path), "status": "present", "sha256": digest,
        "size_bytes": after.st_size, "mtime_ns": after.st_mtime_ns,
    }


def environment_versions() -> dict[str, str]:
    """Record runtime versions relevant to deterministic calculations."""
    # Both platform.platform() and platform.machine() can probe WMI on Windows.
    # Restricted workers may lack WMI and PROCESSOR_* variables; concurrent COM
    # probes can also terminate the interpreter. Use the OS version and Python's
    # build architecture, which do not depend on the worker's account/environment.
    if sys.platform == "win32":
        native = sys.getwindowsversion()
        build_platform = sysconfig.get_platform()
        architecture = {
            "win-amd64": "AMD64", "win32": "x86",
            "win-arm64": "ARM64", "win-arm32": "ARM",
        }.get(build_platform, build_platform)
        os_identity = f"Windows-{native.major}.{native.minor}.{native.build}-{architecture}"
    else:
        os_identity = platform.platform()
    versions = {"python": platform.python_version(), "platform": os_identity}
    for package in ("rdkit", "numpy", "rdchiral", "pypdf"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed"
    return versions


def capture_baseline(
    repository: Path, artifacts: Mapping[str, Path],
    *, include_guidance: bool = True,
    on_progress: Callable[[dict[str, Any]], None] | None = None,
    fingerprint_cache: Path | None = None,
) -> dict[str, Any]:
    """Record a snapshot with optional preparation telemetry and cancellation.

    Callback exceptions abort capture. Telemetry is never part of the baseline
    identity. A configured cache reuses hashes only for unchanged file snapshots;
    omitting it performs full content hashing. Replay always bypasses the cache.
    """
    def report(stage: str, **details: Any) -> None:
        if on_progress is not None:
            on_progress({"stage": stage, **details})

    report("resolving")
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
    expected_hashes: dict[Path, str] = {}
    if "processed_dataset" in paths:
        from condition_recommender.processed_release import resolve_processed_release
        release = resolve_processed_release(paths["processed_dataset"])
        paths["processed_dataset"] = release.root / "manifest.json"
        for name in ("condition_index", "shared_core_index", "fragment_index", "retro_library",
                     "forward_library", "evidence_catalog", "route_catalog", "composite_library", "composite_catalog"):
            if name in release.manifest["artifacts"]:
                selected = release.artifact(name)
                if name in paths and paths[name].resolve() != selected:
                    raise ValueError(f"Artifact {name} conflicts with the selected processed release")
                paths[name] = selected
                expected_hashes[selected] = release.manifest["artifacts"][name]["sha256"]
        for alias in ("procedure_catalog", "reference_catalog"):
            selected = release.artifact("evidence_catalog")
            if alias in paths and paths[alias].resolve() != selected:
                raise ValueError(f"Artifact {alias} conflicts with the selected processed release")
            paths[alias] = selected
            expected_hashes[selected] = release.manifest["artifacts"]["evidence_catalog"]["sha256"]
    if "weak_label_records" in paths:
        from condition_recommender import weak_label_recipe_catalog_path

        catalog = weak_label_recipe_catalog_path(paths["weak_label_records"]).resolve()
        configured = paths.get("weak_label_recipe_catalog", catalog).resolve()
        if configured != catalog:
            raise ValueError("weak_label_recipe_catalog must be the catalog beside weak_label_records")
        paths["weak_label_recipe_catalog"] = catalog
    fingerprints: dict[Path, dict[str, Any]] = {}
    cache = FingerprintCache(fingerprint_cache) if fingerprint_cache is not None else None
    identities = {}
    file_count = len({path.expanduser().resolve() for path in paths.values()})
    for name, path in paths.items():
        resolved = path.expanduser().resolve()
        if resolved not in fingerprints:
            size = resolved.stat().st_size if resolved.is_file() else 0

            def file_progress(bytes_read: int) -> None:
                report("fingerprinting", artifact=name, bytes_read=bytes_read,
                       size_bytes=size, file_index=len(fingerprints) + 1, file_count=file_count)

            file_progress(0)
            fingerprints[resolved] = artifact_identity(
                resolved, on_progress=file_progress if on_progress is not None else None,
                cache=cache, expected_sha256=expected_hashes.get(resolved),
                on_reuse=lambda: report("fingerprint_reused", artifact=name,
                                        file_index=len(fingerprints) + 1, file_count=file_count),
            )
        identities[name] = dict(fingerprints[resolved])
    report("metadata")
    for value in identities.values():
        path = Path(value["path"])
        if value["status"] == "present" and path.suffix == ".sqlite":
            with closing(sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)) as connection:
                tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                if "metadata" in tables:
                    columns = {row[1] for row in connection.execute("PRAGMA table_info(metadata)")}
                    if "payload" in columns:
                        value["metadata"] = [json.loads(row[0]) for row in connection.execute("SELECT payload FROM metadata")]
    report("validating")
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
    from ..agent_context.learning import build_learning_context, disabled_learning_context

    report("guidance")
    baseline["learning_context"] = (build_learning_context(baseline) if include_guidance
                                    else disabled_learning_context())
    baseline["scientific_identity"] = scientific_identity(baseline)
    report("complete")
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
    current_environment = environment_versions()
    recorded_environment = baseline["environment"]
    if current_environment != recorded_environment:
        differences = "; ".join(
            f"{key}: recorded={recorded_environment.get(key)!r}, "
            f"current={current_environment.get(key)!r}"
            for key in sorted(recorded_environment.keys() | current_environment.keys())
            if recorded_environment.get(key) != current_environment.get(key)
        )
        raise ValueError(
            "Scientific runtime changed; align the worker environment and start "
            f"a new investigation ({differences})"
        )
    hashes: dict[Path, str] = {}
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
        if full_hash:
            if path not in hashes:
                hashes[path] = sha256_file(path)
            if hashes[path] != value["sha256"]:
                raise ValueError(f"Baseline artifact content changed: {path}")
