"""Immutable per-turn application context, separate from scientific execution."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .baseline import (
    APPLICATION_DIRECTORY_LAYERS, APPLICATION_MODULE_LAYERS, WORKSPACE_PREFIX,
    scientific_identity, sha256_file, verify_baseline,
)
from .store import InvestigationEvent, InvestigationStore, canonical_bytes


APPLICATION_CONTEXT_SCHEMA = "scientific_application_context.v1"
PRESENTATION_FILES = (
    "app/web_api/scientific_presentation.py",
    "app/web_api/scientific_chat.js",
    "app/web_api/scientific_chat.css",
    "app/web_api/scientific_chat.html",
)


def _digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _resource_directories(repository: Path) -> list[tuple[str, str, Path]]:
    directory = repository / "chem_coworker" / "scientific_workspace"
    return [
        (name, layer, directory / name if (directory / name).is_dir()
         else Path(__file__).parent / name)
        for name, layer in APPLICATION_DIRECTORY_LAYERS.items()
    ]


def application_manifest(repository: Path) -> dict[str, dict[str, str]]:
    """Identify runtime, guidance, and presentation files by explicit ownership."""
    directory = repository / "chem_coworker" / "scientific_workspace"
    layers: dict[str, dict[str, str]] = {
        "runtime": {}, "guidance": {}, "presentation": {},
    }
    for name, layer in APPLICATION_MODULE_LAYERS.items():
        path = directory / name
        if path.is_file():
            layers[layer][WORKSPACE_PREFIX + name] = sha256_file(path)
    for name, layer, source in _resource_directories(repository):
        for path in sorted(source.rglob("*.md")):
            key = WORKSPACE_PREFIX + name + "/" + path.relative_to(source).as_posix()
            layers[layer][key] = sha256_file(path)
    for relative in PRESENTATION_FILES:
        path = repository / relative
        if path.is_file():
            layers["presentation"][relative] = sha256_file(path)
    return layers


def _resources(repository: Path, layers: Mapping[str, Mapping[str, str]]) -> dict[str, Any]:
    resources = {}
    for name, layer, source in _resource_directories(repository):
        for path in sorted(source.rglob("*.md")):
            content = path.read_bytes()
            digest = hashlib.sha256(content).hexdigest()
            key = name + "/" + path.relative_to(source).as_posix()
            if digest != layers[layer].get(WORKSPACE_PREFIX + key):
                raise ValueError(f"Application resource changed during capture: {path}")
            resources[key] = {"sha256": digest, "text": content.decode("utf-8")}
    expected = {
        key[len(WORKSPACE_PREFIX):]
        for files in layers.values() for key in files
        if key.startswith(WORKSPACE_PREFIX) and key.endswith(".md")
    }
    if expected != set(resources):
        raise ValueError("Application resource inventory changed during capture")
    return resources


def current_application_context(store: InvestigationStore) -> dict[str, Any] | None:
    """Read the latest recorded turn context without accepting unrecorded changes."""
    for event in reversed(store.events()):
        if event.kind != "application_context":
            continue
        context = store.read_artifact(event.artifact_ref)
        if context.get("schema_version") != APPLICATION_CONTEXT_SCHEMA:
            raise ValueError("Unsupported application context schema")
        payload = {key: value for key, value in context.items() if key != "sha256"}
        if _digest(payload) != context.get("sha256"):
            raise ValueError("Application context checksum mismatch")
        return context
    return None


def capture_application_context(
    store: InvestigationStore, *, runtime: Mapping[str, Any],
) -> dict[str, Any]:
    """Snapshot application resources and pinned advice for the next agent turn.

    A changed task guide or runtime configuration may enter only through an
    explicit context record. Scientific changes still require a new baseline.
    The final composed prompt is saved separately by conversation orchestration.
    """
    from .learning import build_learning_context

    baseline = store.manifest["baseline"]
    verify_baseline(baseline)
    repository = Path(baseline["repository"])
    layers = application_manifest(repository)
    resources = _resources(repository, layers)
    previous = current_application_context(store)
    pinned = (previous or {}).get("learning_context") or baseline.get("learning_context")
    guides = {
        Path(key).stem: {"task": Path(key).stem, **resource}
        for key, resource in resources.items()
        if key.startswith("guides/") and len(Path(key).parts) == 2
    }
    if pinned is None:
        learning = build_learning_context(
            baseline,
            guide_snapshots=(guides if (repository / WORKSPACE_PREFIX / "guides").is_dir()
                             else None),
        )
    else:
        # Advice remains pinned for the investigation. Updating a guide must not
        # silently admit lessons published during an earlier turn of this run.
        learning_payload = {key: value for key, value in pinned.items() if key != "sha256"}
        if _digest(learning_payload) != pinned.get("sha256"):
            raise ValueError("Learning context checksum mismatch")
        if (repository / WORKSPACE_PREFIX / "guides").is_dir():
            learning_payload["guides"] = guides
        learning = {**learning_payload, "sha256": _digest(learning_payload)}
    # Parse a canonical snapshot so callers cannot mutate nested runtime settings
    # after capture. Only JSON data is recorded; it is never loaded as code.
    runtime_snapshot = json.loads(canonical_bytes(dict(runtime)))
    payload = {
        "schema_version": APPLICATION_CONTEXT_SCHEMA,
        "scientific_identity": scientific_identity(baseline),
        "layers": {
            layer: {"code_files": files, "sha256": _digest(files)}
            for layer, files in layers.items()
        },
        "resources": resources,
        "resource_directories": {
            name: str(source) for name, _, source in _resource_directories(repository)
        },
        "learning_context": learning,
        "runtime": runtime_snapshot,
        "guidance_identity": "sha256:" + _digest({
            "guidance_files": layers["guidance"],
            "presentation_files": layers["presentation"],
            "resources": resources, "learning_context": learning,
        }),
        "authority": "application_configuration_not_scientific_evidence",
    }
    verify_baseline(baseline)
    # The pinned baseline may supply nested lesson lists. Return an independent
    # JSON snapshot so modifying a capture cannot mutate the baseline in memory.
    snapshot = json.loads(canonical_bytes(payload))
    return {**snapshot, "sha256": _digest(snapshot)}


def record_application_context(
    store: InvestigationStore, *, turn_id: str, runtime: Mapping[str, Any],
) -> InvestigationEvent:
    """Record the context used by one turn without modifying its investigation."""
    if not isinstance(turn_id, str) or not turn_id.strip():
        raise ValueError("An application context requires a turn_id")
    context = capture_application_context(store, runtime=runtime)
    payload = {key: value for key, value in context.items() if key != "sha256"}
    payload["turn_id"] = turn_id
    return store.append("application_context", {**payload, "sha256": _digest(payload)})
