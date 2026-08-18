"""Versioned admission policy for executable generic reaction cores."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal, Mapping


CoreAdmissionPolicyName = Literal["pass_only", "validated_departures"]
CORE_ADMISSION_POLICY_NAMES: tuple[CoreAdmissionPolicyName, ...] = (
    "pass_only",
    "validated_departures",
)
_DEFINITION_PATH = (
    Path(__file__).resolve().parent
    / "definitions"
    / "generic_core_admission_policy.v1.json"
)


@dataclass(frozen=True)
class GenericCoreAdmissionPolicy:
    """Declarative threshold for sending a reaction core to compilation."""

    name: CoreAdmissionPolicyName
    definition_id: str
    allowed_core_statuses: tuple[str, ...]
    allowed_review_reasons: tuple[str, ...]
    allowed_edit_evidence: tuple[str, ...]
    minimum_active_atom_mapping_coverage: float
    required_completeness_status: str
    required_passed_checks: tuple[str, ...]
    require_validated_departures_for_review: bool

    def __post_init__(self) -> None:
        if self.name not in CORE_ADMISSION_POLICY_NAMES:
            raise ValueError("unsupported generic core-admission policy")
        if not self.definition_id.startswith(
            "generic_core_admission_policy.v1@"
        ):
            raise ValueError("unexpected generic core-admission definition")
        if not self.allowed_core_statuses or "blocked" in self.allowed_core_statuses:
            raise ValueError("blocked cores cannot be admitted")
        if not 0.0 <= self.minimum_active_atom_mapping_coverage <= 1.0:
            raise ValueError("core-admission mapping threshold is invalid")
        if self.required_completeness_status != "verified":
            raise ValueError("executable cores require verified completeness")
        if self.name == "pass_only" and self.allowed_core_statuses != ("pass",):
            raise ValueError("pass-only policy may admit only pass cores")
        if self.require_validated_departures_for_review and not all(
            (
                "review" in self.allowed_core_statuses,
                self.allowed_review_reasons,
                self.allowed_edit_evidence,
                self.required_passed_checks,
            )
        ):
            raise ValueError("departure review policy lacks evidence constraints")


@dataclass(frozen=True)
class GenericCoreAdmissionDecision:
    """One transparent compiler core-admission decision."""

    admitted: bool
    policy_name: CoreAdmissionPolicyName
    policy_definition_id: str
    reason: str
    checked_edit_count: int = 0
    departure_edit_descriptors: tuple[str, ...] = ()
    issues: tuple[str, ...] = ()

    def diagnostics(self) -> dict[str, Any]:
        """Return stable JSON-compatible compiler diagnostics."""

        value = asdict(self)
        value["departure_edit_descriptors"] = list(
            self.departure_edit_descriptors
        )
        value["issues"] = list(self.issues)
        return {f"core_admission_{key}": item for key, item in value.items()}


@lru_cache(maxsize=1)
def _load_policy_values() -> dict[str, Any]:
    value = json.loads(_DEFINITION_PATH.read_text(encoding="utf-8"))
    if value.get("schema_version") != "1.0":
        raise ValueError("unsupported generic core-admission schema")
    if value.get("definition_id") != "generic_core_admission_policy.v1":
        raise ValueError("unexpected generic core-admission definition ID")
    if not str(value.get("definition_version") or ""):
        raise ValueError("generic core-admission definition requires a version")
    if set(value.get("policies") or {}) != set(CORE_ADMISSION_POLICY_NAMES):
        raise ValueError("generic core-admission policies are incomplete")
    return value


@lru_cache(maxsize=len(CORE_ADMISSION_POLICY_NAMES))
def load_generic_core_admission_policy(
    name: CoreAdmissionPolicyName = "pass_only",
) -> GenericCoreAdmissionPolicy:
    """Load and validate one versioned generic core-admission policy."""

    if name not in CORE_ADMISSION_POLICY_NAMES:
        raise ValueError("unsupported generic core-admission policy")
    value = _load_policy_values()
    raw = dict(value["policies"][name])
    return GenericCoreAdmissionPolicy(
        name=name,
        definition_id=(
            f"{value['definition_id']}@{value['definition_version']}:{name}"
        ),
        allowed_core_statuses=tuple(
            str(item) for item in raw.get("allowed_core_statuses") or ()
        ),
        allowed_review_reasons=tuple(
            str(item) for item in raw.get("allowed_review_reasons") or ()
        ),
        allowed_edit_evidence=tuple(
            str(item) for item in raw.get("allowed_edit_evidence") or ()
        ),
        minimum_active_atom_mapping_coverage=float(
            raw["minimum_active_atom_mapping_coverage"]
        ),
        required_completeness_status=str(
            raw["required_completeness_status"]
        ),
        required_passed_checks=tuple(
            str(item) for item in raw.get("required_passed_checks") or ()
        ),
        require_validated_departures_for_review=bool(
            raw["require_validated_departures_for_review"]
        ),
    )


def _side_map_numbers(
    observation: Mapping[str, Any],
    field: str,
) -> set[int]:
    values: set[int] = set()
    for participant in observation.get(field) or ():
        structure = (participant or {}).get("molecular_structure") or {}
        for component in structure.get("components") or ():
            for atom in (component or {}).get("atoms") or ():
                map_number = int((atom or {}).get("atom_map_number") or 0)
                if map_number > 0:
                    values.add(map_number)
    return values


def _edit_map_number(edit: Mapping[str, Any], field: str) -> int:
    return int(((edit.get(field) or {}).get("atom_map_number") or 0))


def _departure_descriptor(
    edit: Mapping[str, Any],
    departing_field: str,
) -> str:
    retained_field = "atom_2" if departing_field == "atom_1" else "atom_1"
    retained = edit.get(retained_field) or {}
    departing = edit.get(departing_field) or {}
    return ":".join(
        (
            str(edit.get("edit_type") or "unknown"),
            f"{retained.get('element') or '?'}-{departing.get('element') or '?'}",
            str(edit.get("old_order") or "NONE"),
            str(_edit_map_number(edit, retained_field)),
            str(_edit_map_number(edit, departing_field)),
        )
    )


def _review_departure_evidence(
    core: Mapping[str, Any],
    observation: Mapping[str, Any],
    policy: GenericCoreAdmissionPolicy,
) -> tuple[int, tuple[str, ...], tuple[str, ...]]:
    quality = core.get("quality") or {}
    edits = tuple(observation.get("edits") or ())
    reactant_maps = _side_map_numbers(observation, "reactants")
    product_maps = _side_map_numbers(observation, "products")
    issues: list[str] = []
    departures: list[str] = []
    checked = 0
    if len(edits) != int(quality.get("edit_count") or 0):
        issues.append("edit_count_mismatch")
    heavy_count = sum((edit or {}).get("atom_2") is not None for edit in edits)
    if heavy_count != int(quality.get("heavy_atom_edit_count") or 0):
        issues.append("heavy_atom_edit_count_mismatch")
    if not reactant_maps or not product_maps:
        issues.append("mapped_participant_graphs_unavailable")
    for raw_edit in edits:
        edit = raw_edit or {}
        first = _edit_map_number(edit, "atom_1")
        second_present = edit.get("atom_2") is not None
        second = _edit_map_number(edit, "atom_2") if second_present else 0
        evidence = str(edit.get("evidence") or "")
        confidence = float(edit.get("confidence") or 0.0)
        if evidence not in policy.allowed_edit_evidence or confidence < 1.0:
            issues.append("edit_mapping_evidence_not_validated")
            continue
        if first <= 0 or first not in reactant_maps:
            issues.append("edit_atom_missing_from_mapped_reactants")
            continue
        if not second_present:
            if first in product_maps:
                checked += 1
            else:
                issues.append("hydrogen_carrier_missing_from_product")
            continue
        if second <= 0 or second not in reactant_maps:
            issues.append("edit_atom_missing_from_mapped_reactants")
            continue
        first_retained = first in product_maps
        second_retained = second in product_maps
        if first_retained and second_retained:
            checked += 1
            continue
        is_departure = bool(
            str(edit.get("edit_type") or "") == "broken"
            and edit.get("old_order") is not None
            and edit.get("new_order") is None
            and first_retained != second_retained
        )
        if not is_departure:
            issues.append("unchecked_edit_is_not_validated_departure")
            continue
        departing_field = "atom_2" if first_retained else "atom_1"
        departures.append(_departure_descriptor(edit, departing_field))
    reported_checked = float(quality.get("checked_edit_fraction") or 0.0) * int(
        quality.get("edit_count") or 0
    )
    if not math.isclose(reported_checked, checked, abs_tol=1e-9):
        issues.append("checked_edit_count_mismatch")
    if not departures:
        issues.append("validated_departure_missing")
    return checked, tuple(sorted(set(departures))), tuple(sorted(set(issues)))


def assess_generic_core_admission(
    core: Mapping[str, Any],
    observation: Mapping[str, Any],
    *,
    policy_name: CoreAdmissionPolicyName = "pass_only",
) -> GenericCoreAdmissionDecision:
    """Decide whether one core may proceed to executable compilation."""

    policy = load_generic_core_admission_policy(policy_name)
    quality = core.get("quality") or {}
    status = str(quality.get("status") or "missing")
    if status == "pass":
        return GenericCoreAdmissionDecision(
            admitted=True,
            policy_name=policy.name,
            policy_definition_id=policy.definition_id,
            reason="pass_core",
            checked_edit_count=int(quality.get("edit_count") or 0),
        )
    if status not in policy.allowed_core_statuses:
        return GenericCoreAdmissionDecision(
            admitted=False,
            policy_name=policy.name,
            policy_definition_id=policy.definition_id,
            reason="core_status_not_allowed",
        )
    review_reasons = set(quality.get("review_reasons") or ())
    required_checks = set(policy.required_passed_checks)
    issues = []
    if not review_reasons or not review_reasons.issubset(
        policy.allowed_review_reasons
    ):
        issues.append("unsupported_core_review_reason")
    if quality.get("blocking_reasons"):
        issues.append("core_has_blocking_reasons")
    if float(quality.get("active_atom_mapping_coverage") or 0.0) < (
        policy.minimum_active_atom_mapping_coverage
    ):
        issues.append("active_atom_mapping_below_policy")
    if not required_checks.issubset(quality.get("passed_checks") or ()):
        issues.append("required_structural_checks_missing")
    completeness = observation.get("completeness") or {}
    if completeness.get("status") != policy.required_completeness_status:
        issues.append("product_completeness_not_verified")
    checked = 0
    departures: tuple[str, ...] = ()
    if policy.require_validated_departures_for_review:
        checked, departures, departure_issues = _review_departure_evidence(
            core,
            observation,
            policy,
        )
        issues.extend(departure_issues)
    unique_issues = tuple(sorted(set(issues)))
    return GenericCoreAdmissionDecision(
        admitted=not unique_issues,
        policy_name=policy.name,
        policy_definition_id=policy.definition_id,
        reason=(
            "validated_departure_review"
            if not unique_issues
            else "review_core_policy_rejected"
        ),
        checked_edit_count=checked,
        departure_edit_descriptors=departures,
        issues=unique_issues,
    )


__all__ = [
    "CORE_ADMISSION_POLICY_NAMES",
    "CoreAdmissionPolicyName",
    "GenericCoreAdmissionDecision",
    "GenericCoreAdmissionPolicy",
    "assess_generic_core_admission",
    "load_generic_core_admission_policy",
]
