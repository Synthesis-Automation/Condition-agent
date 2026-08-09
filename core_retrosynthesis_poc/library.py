"""Build and persist core-derived retrosynthesis template libraries."""

from __future__ import annotations

import gzip
import json
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, TextIO

from retrosynthesis_poc.library import iter_rows

from .compiler import compile_core_templates
from .models import (
    AbstractionLevel,
    CoreLibraryBuildReport,
    CoreTemplate,
    CoreTemplateLibrary,
    CoreTemplatePrecedent,
)


LIBRARY_DEFINITION = {
    "definition_id": "core_retrosynthesis_poc.v1",
    "allowed_bonds": ["C-N", "C-O", "C-S"],
    "abstraction_levels": ["L1", "L2"],
    "L1": "edit atoms plus connected precursor-only handle subgraphs",
    "L2": "L1 chemistry plus the first molecular shell",
    "source_round_trip_required": True,
    "context_embedded_in_smarts": False,
    "context_features": [
        "reaction_core_substituent_profiles",
        "spectator_groups",
        "steric_accessibility",
        "electronic_activation",
    ],
}


@dataclass
class _Accumulator:
    template: CoreTemplate
    observation_support: int = 0
    support_units: set[str] = field(default_factory=set)
    precedents: list[CoreTemplatePrecedent] = field(default_factory=list)


def _open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def build_library(
    rows: Iterable[Dict[str, Any]],
    *,
    levels: Iterable[AbstractionLevel] = ("L1", "L2"),
    max_rows: Optional[int] = None,
    max_precedents_per_template: int = 8,
) -> tuple[CoreTemplateLibrary, CoreLibraryBuildReport]:
    """Compile, round-trip, deduplicate, and aggregate core-derived rules."""

    if max_rows is not None and max_rows < 1:
        raise ValueError("max_rows must be positive")
    if max_precedents_per_template < 1:
        raise ValueError("max precedents must be positive")
    selected_levels = tuple(dict.fromkeys(levels))
    accumulators: dict[str, _Accumulator] = {}
    rejections: Counter[str] = Counter()
    source_count = 0
    accepted_count = 0
    for row in rows:
        if max_rows is not None and source_count >= max_rows:
            break
        source_count += 1
        result = compile_core_templates(row, levels=selected_levels)
        if not result.templates:
            rejections[str(result.rejection_reason or "unknown_rejection")] += 1
            continue
        accepted_count += 1
        for template in result.templates:
            accumulator = accumulators.setdefault(
                template.template_id,
                _Accumulator(template=template),
            )
            accumulator.observation_support += 1
            precedent = template.precedents[0]
            accumulator.support_units.add(precedent.support_unit_id)
            if (
                len(accumulator.precedents) < max_precedents_per_template
                and precedent not in accumulator.precedents
            ):
                accumulator.precedents.append(precedent)

    provisional = []
    for template_id, accumulator in sorted(accumulators.items()):
        del template_id
        provisional.append(
            replace(
                accumulator.template,
                observation_support=accumulator.observation_support,
                independent_reference_support=len(accumulator.support_units),
                precedents=tuple(accumulator.precedents),
            )
        )
    operator_observations: Counter[str] = Counter()
    operator_references: dict[str, set[str]] = {}
    for accumulator in accumulators.values():
        operator_id = accumulator.template.operator_id
        operator_observations[operator_id] += accumulator.observation_support
        operator_references.setdefault(operator_id, set()).update(
            accumulator.support_units
        )
    templates = tuple(
        replace(
            template,
            operator_observation_support=operator_observations[template.operator_id],
            operator_reference_support=len(operator_references[template.operator_id]),
        )
        for template in sorted(provisional, key=lambda value: value.template_id)
    )
    rejection_counts = dict(sorted(rejections.items()))
    library = CoreTemplateLibrary(
        templates=templates,
        source_row_count=source_count,
        accepted_observation_count=accepted_count,
        rejection_counts=rejection_counts,
        definition={**LIBRARY_DEFINITION, "compiled_levels": list(selected_levels)},
    )
    report = CoreLibraryBuildReport(
        source_row_count=source_count,
        accepted_observation_count=accepted_count,
        unique_template_count=len(templates),
        unique_operator_count=len({template.operator_id for template in templates}),
        rejection_counts=rejection_counts,
    )
    return library, report


def save_library(library: CoreTemplateLibrary, destination: str | Path) -> None:
    """Write canonical JSON, optionally gzip-compressed by filename suffix."""

    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        library.to_dict(),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
            handle.write(payload)
            handle.write("\n")
    else:
        path.write_text(payload + "\n", encoding="utf-8")


def load_library(source: str | Path) -> CoreTemplateLibrary:
    """Load and validate a core-derived template library."""

    path = Path(source)
    with _open_text(path) as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("core retrosynthesis library must be a JSON object")
    return CoreTemplateLibrary.from_dict(value)


__all__ = [
    "LIBRARY_DEFINITION",
    "build_library",
    "iter_rows",
    "load_library",
    "save_library",
]
