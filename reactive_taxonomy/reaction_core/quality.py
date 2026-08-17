"""Deterministic validation and quality assessment for reaction cores."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Sequence, Tuple

from ..reaction_models import ReactionEdit
from ..reaction_edits import PROJECTED_UNMAPPED_DEPARTING_BOUNDARY_EVIDENCE
from .common import Location
from .models import ReactionCoreQuality


_RULES_PATH = (
    Path(__file__).parents[1]
    / "definitions"
    / "reaction_core_quality.v1.json"
)


@lru_cache(maxsize=1)
def load_reaction_core_quality_rules() -> dict[str, Any]:
    """Load and validate the versioned reaction-core quality policy."""
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        rules = dict(json.load(handle))
    if str(rules.get("schema_version") or "") != "1.0":
        raise ValueError("unsupported reaction-core quality schema")
    if str(rules.get("definition_id") or "") != "reaction_core_quality.v1":
        raise ValueError("unexpected reaction-core quality definition ID")
    if not str(rules.get("definition_version") or ""):
        raise ValueError("reaction-core quality definition requires a version")
    coverage = float(rules["minimum_active_atom_mapping_coverage"])
    if not 0.0 <= coverage <= 1.0:
        raise ValueError("reaction-core mapping coverage threshold is invalid")
    if int(rules["maximum_active_atom_count_without_review"]) < 1:
        raise ValueError("reaction-core active atom threshold must be positive")
    if int(rules["maximum_event_count_without_review"]) < 1:
        raise ValueError("reaction-core event threshold must be positive")
    return rules


def _mapped_location(
    edit_atom: Any,
    locations: Mapping[int, Location],
) -> Location | None:
    map_number = edit_atom.atom_map_number
    return locations.get(int(map_number)) if map_number is not None else None


def _bond_order(
    first: Location,
    second: Location,
) -> str | None:
    first_component, first_molecule, first_index = first
    second_component, _, second_index = second
    if first_component.component_index != second_component.component_index:
        return None
    bond = first_molecule.GetBondBetweenAtoms(first_index, second_index)
    return str(bond.GetBondType()).upper() if bond is not None else None


def validate_core_edits(
    edits: Sequence[ReactionEdit],
    *,
    reactant_by_map: Mapping[int, Location],
    product_by_map: Mapping[int, Location],
) -> tuple[int, Tuple[str, ...]]:
    """Check mapped edits against the supplied reactant and product graphs."""
    checked = 0
    issues = []
    for edit in edits:
        reactant_1 = _mapped_location(edit.atom_1, reactant_by_map)
        product_1 = _mapped_location(edit.atom_1, product_by_map)
        if edit.atom_2 is None:
            if reactant_1 is None or product_1 is None:
                continue
            checked += 1
            before_h = int(
                reactant_1[1]
                .GetAtomWithIdx(reactant_1[2])
                .GetTotalNumHs(includeNeighbors=True)
            )
            after_h = int(
                product_1[1]
                .GetAtomWithIdx(product_1[2])
                .GetTotalNumHs(includeNeighbors=True)
            )
            if edit.old_order is not None and after_h >= before_h:
                issues.append("hydrogen_loss_inconsistent")
            if edit.new_order is not None and after_h <= before_h:
                issues.append("hydrogen_gain_inconsistent")
            continue

        reactant_2 = _mapped_location(edit.atom_2, reactant_by_map)
        product_2 = _mapped_location(edit.atom_2, product_by_map)
        if any(
            value is None
            for value in (reactant_1, reactant_2, product_1, product_2)
        ):
            if (
                edit.evidence
                == PROJECTED_UNMAPPED_DEPARTING_BOUNDARY_EVIDENCE
                and edit.edit_type == "broken"
                and edit.old_order is not None
                and edit.new_order is None
                and reactant_1 is not None
                and product_1 is not None
                and edit.atom_2.atom_map_number is None
            ):
                checked += 1
            continue
        checked += 1
        before_order = _bond_order(reactant_1, reactant_2)  # type: ignore[arg-type]
        after_order = _bond_order(product_1, product_2)  # type: ignore[arg-type]
        expected_before = (
            str(edit.old_order).upper() if edit.old_order is not None else None
        )
        expected_after = (
            str(edit.new_order).upper() if edit.new_order is not None else None
        )
        if before_order != expected_before or after_order != expected_after:
            issues.append(f"{edit.edit_type}_bond_state_inconsistent")
    return checked, tuple(sorted(set(issues)))


def assess_reaction_core_quality(
    *,
    active_atom_count: int,
    mapped_active_atom_count: int,
    edit_count: int,
    heavy_atom_edit_count: int,
    checked_edit_count: int,
    consistency_issues: Sequence[str],
    event_count: int,
    remote_continuity_unresolved: bool,
    no_op_primary_center: bool,
    unmapped_active_atoms_are_validated_departures: bool = False,
) -> ReactionCoreQuality:
    """Build a transparent pass/review/blocked core-quality assessment."""
    rules = load_reaction_core_quality_rules()
    coverage = (
        mapped_active_atom_count / active_atom_count
        if active_atom_count
        else 0.0
    )
    checked_fraction = checked_edit_count / edit_count if edit_count else 0.0
    passed = []
    review = []
    blocked = list(consistency_issues)

    if no_op_primary_center:
        blocked.append("no_op_primary_center")
    else:
        passed.append("primary_center_changes_state")
    if event_count < 1:
        blocked.append("reaction_core_event_missing")
    elif event_count > int(rules["maximum_event_count_without_review"]):
        review.append("many_disconnected_edit_events")
    else:
        passed.append("event_count_within_review_limit")
    if active_atom_count > int(
        rules["maximum_active_atom_count_without_review"]
    ):
        review.append("large_active_core")
    else:
        passed.append("active_core_size_within_review_limit")
    if (
        coverage < float(rules["minimum_active_atom_mapping_coverage"])
        and not unmapped_active_atoms_are_validated_departures
    ):
        review.append("partial_active_atom_mapping")
    elif coverage < float(rules["minimum_active_atom_mapping_coverage"]):
        passed.append("partial_mapping_limited_to_validated_departures")
    else:
        passed.append("active_atom_mapping_complete")
    if checked_fraction < 1.0:
        review.append("not_all_edits_graph_checked")
    elif not consistency_issues:
        passed.append("all_edits_graph_consistent")
    if remote_continuity_unresolved:
        review.append("remote_continuity_unresolved")
    else:
        passed.append("remote_continuity_resolved")

    blocked_values = tuple(sorted(set(blocked)))
    review_values = tuple(sorted(set(review)))
    status = (
        "blocked"
        if blocked_values
        else "review"
        if review_values
        else "pass"
    )
    return ReactionCoreQuality(
        status=status,
        active_atom_mapping_coverage=coverage,
        checked_edit_fraction=checked_fraction,
        edit_count=edit_count,
        heavy_atom_edit_count=heavy_atom_edit_count,
        event_count=event_count,
        passed_checks=tuple(sorted(set(passed))),
        review_reasons=review_values,
        blocking_reasons=blocked_values,
        definition_version=(
            f"{rules['definition_id']}@{rules['definition_version']}"
        ),
    )


__all__ = [
    "assess_reaction_core_quality",
    "load_reaction_core_quality_rules",
    "validate_core_edits",
]
