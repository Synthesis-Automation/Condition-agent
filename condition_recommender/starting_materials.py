"""Evidence-preserving terminal-material assessment for agent-owned planning."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from math import isfinite
from pathlib import Path
import sqlite3
from typing import Any

from condition_registry.api import get_registry
from condition_registry.resolver import ConditionRegistry
from condition_registry.vocabulary import condition_registry_definition_versions
from reactive_taxonomy.material_identity import MaterialIdentity, identify_material

from .fragment_search import lookup_exact_product


POLICY_PATH = Path(__file__).with_name("definitions") / "starting_material_policy.v1.json"


def _positive_number(value: Any, name: str) -> float:
    if type(value) not in (int, float) or not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return float(value)


@dataclass(frozen=True)
class StartingMaterialPolicy:
    """Validated, versioned stopping rules; definitions never name executable code."""

    schema_version: str
    definition_version: str
    default_mw_threshold: float
    precedent_limit: int
    check_order: tuple[str, ...]
    mw_comparison: str
    availability: str

    def __post_init__(self) -> None:
        if self.schema_version != "starting_material_policy.v1":
            raise ValueError("Unsupported starting-material policy schema")
        if not isinstance(self.definition_version, str) or not self.definition_version.strip():
            raise ValueError("Starting-material policy requires a definition version")
        _positive_number(self.default_mw_threshold, "default_mw_threshold")
        if type(self.precedent_limit) is not int or not 1 <= self.precedent_limit <= 10:
            raise ValueError("precedent_limit must be an integer between 1 and 10")
        if self.check_order != ("registry", "literature", "molecular_weight"):
            raise ValueError("Unsupported starting-material check order")
        if self.mw_comparison != "strictly_less_than" or self.availability != "unknown":
            raise ValueError("Unsupported stopping comparison or availability claim")


def load_starting_material_policy(path: str | Path = POLICY_PATH) -> StartingMaterialPolicy:
    """Load and validate the complete policy definition, rejecting unknown fields."""
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    expected = set(StartingMaterialPolicy.__dataclass_fields__)
    if not isinstance(value, dict) or set(value) != expected:
        raise ValueError("Invalid starting-material policy fields")
    if not isinstance(value["check_order"], list):
        raise ValueError("check_order must be a list")
    value["check_order"] = tuple(value["check_order"])
    return StartingMaterialPolicy(**value)


@dataclass(frozen=True)
class StartingMaterialAssessment:
    """A planning decision and its assumptions, never confirmation of supply."""

    input_smiles: str
    canonical_smiles: str
    material: MaterialIdentity
    molecular_weight: float
    mw_threshold: float
    stop_expansion: bool
    stop_reason: str
    status: str
    availability: str
    registry: dict[str, Any]
    literature: dict[str, Any]
    policy: dict[str, Any]
    definition_version: str
    unavailable_starting_materials: tuple[str, ...]
    warnings: tuple[str, ...]
    schema_version: str = "starting_material_assessment.v1"

    def to_dict(self) -> dict[str, Any]:
        """Serialize all checks, including skipped and failed checks."""
        return asdict(self)


def _registry_lookup(canonical: str, registry: ConditionRegistry | None) -> dict[str, Any]:
    try:
        selected = registry if registry is not None else get_registry()
        resolved = selected.resolve_identifier(canonical, identifier_type="smiles")
        evidence = resolved.to_dict()
        evidence["source"] = "condition_registry" if registry is None else "supplied_registry"
        evidence["definition_versions"] = (
            dict(condition_registry_definition_versions()) if registry is None else {}
        )
        if resolved.substance is not None and resolved.substance.status != "active":
            evidence["status"] = "inactive"
        return evidence
    except (OSError, ValueError) as exc:
        return {"status": "error", "error_type": type(exc).__name__, "message": str(exc)}


def _literature_lookup(canonical: str, path: str | Path | None, limit: int) -> dict[str, Any]:
    if path is None:
        return {"status": "unavailable", "message": "No fragment_index configured"}
    try:
        return lookup_exact_product(path, canonical, limit=limit)
    except FileNotFoundError as exc:
        return {"status": "unavailable", "error_type": type(exc).__name__, "message": str(exc)}
    except (OSError, ValueError, sqlite3.Error) as exc:
        return {"status": "error", "error_type": type(exc).__name__, "message": str(exc)}


def assess_starting_material(
    smiles: str, *, fragment_index: str | Path | None = None,
    mw_threshold: float | None = None, allow_registry_stop: bool = True,
    allow_literature_stop: bool = True, allow_mw_stop: bool = True,
    unavailable_starting_materials: list[str] | tuple[str, ...] = (),
    registry: ConditionRegistry | None = None,
) -> StartingMaterialAssessment:
    """Check registry identity, exact reported products, then a strict MW cutoff.

    Disabled stages are skipped. Explicit exclusions override all stopping rules.
    Ambiguous/inactive registry identities cannot trigger the registry stop.
    Lookup failures remain visible even when a later heuristic permits stopping.
    Every accepted leaf is an assumption requiring availability confirmation.
    """
    policy = load_starting_material_policy()
    threshold = _positive_number(
        policy.default_mw_threshold if mw_threshold is None else mw_threshold, "mw_threshold",
    )
    flags = {"allow_registry_stop": allow_registry_stop,
             "allow_literature_stop": allow_literature_stop, "allow_mw_stop": allow_mw_stop}
    if any(type(value) is not bool for value in flags.values()):
        raise ValueError("Stopping options must be booleans")
    if not isinstance(unavailable_starting_materials, (list, tuple)):
        raise ValueError("unavailable_starting_materials must be a list or tuple of SMILES")
    material = identify_material(smiles)
    excluded = tuple(sorted({identify_material(value).canonical_smiles
                             for value in unavailable_starting_materials}))
    registry_evidence: dict[str, Any] = {"status": "not_run", "reason": "earlier_decision"}
    literature_evidence: dict[str, Any] = {"status": "not_run", "reason": "earlier_decision"}
    warnings = list(material.warnings)
    stop = False
    reason = "no_stopping_criterion"

    if material.canonical_smiles in excluded:
        reason = "explicitly_unavailable"
        warnings.append("Explicitly unavailable material; all stopping rules were bypassed.")
    else:
        registry_evidence = (
            _registry_lookup(material.canonical_smiles, registry) if allow_registry_stop
            else {"status": "not_run", "reason": "disabled_by_policy"}
        )
        status = registry_evidence["status"]
        if status == "resolved":
            stop, reason = True, "registry_match"
            warnings.append("Obtainability is assumed from registry membership; commercial availability is unverified.")
        elif status in {"ambiguous", "inactive", "error"}:
            warnings.append(f"Registry lookup is {status}; no registry-based stopping evidence was accepted.")

        if not stop:
            literature_evidence = (
                _literature_lookup(material.canonical_smiles, fragment_index, policy.precedent_limit)
                if allow_literature_stop else {"status": "not_run", "reason": "disabled_by_policy"}
            )
            status = literature_evidence["status"]
            if status == "matched":
                stop, reason = True, "exact_literature_match"
                warnings.append("Expansion stopped at a reported product; inspect its preparation and confirm obtainability.")
            elif status in {"unavailable", "error"}:
                warnings.append(f"Literature lookup is {status}; absence of this molecule has not been established.")
            warnings.extend(literature_evidence.get("warnings", ()))

        if not stop and allow_mw_stop and material.molecular_weight < threshold:
            stop, reason = True, "low_molecular_weight"
            warnings.append("Expansion stopped by molecular-weight policy; availability and ease of synthesis are unverified.")

    return StartingMaterialAssessment(
        input_smiles=smiles, canonical_smiles=material.canonical_smiles, material=material,
        molecular_weight=material.molecular_weight, mw_threshold=threshold,
        stop_expansion=stop, stop_reason=reason,
        status="assumed_terminal" if stop else "excluded" if reason == "explicitly_unavailable" else "unresolved",
        availability=policy.availability, registry=registry_evidence, literature=literature_evidence,
        policy={**asdict(policy), **flags, "effective_mw_threshold": threshold},
        definition_version=policy.definition_version,
        unavailable_starting_materials=excluded, warnings=tuple(warnings),
    )
