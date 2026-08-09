"""Build and persist structurally diverse generic template libraries."""

from __future__ import annotations

import gzip
import json
from collections import Counter
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, Iterable, Literal

from .generic_compiler import compile_generic_templates
from .generic_models import GenericCoreTemplate, GenericTemplateLibrary


GENERIC_LIBRARY_DEFINITION = {
    "definition_id": "generic_core_retrosynthesis_poc.v1",
    "routing_source": "normalized_graph_edits",
    "named_families_used_for_routing": False,
    "source_round_trip_required": True,
    "supported_transformations": [
        "acyl_substitution",
        "c_c_coupling",
        "carbonyl_oxidation",
        "carbonyl_reduction",
        "conjugate_addition",
        "carbonyl_condensation",
        "ring_formation",
    ],
}


def build_generic_library(
    rows: Iterable[Dict[str, Any]],
    *,
    engine: Literal["reaction_core", "rdchiral"] = "reaction_core",
    levels: Iterable[Literal["L1", "L2"]] = ("L1", "L2"),
    max_precedents_per_template: int = 8,
) -> GenericTemplateLibrary:
    """Compile and aggregate a generic source-round-tripped library."""

    if max_precedents_per_template < 1:
        raise ValueError("max precedents must be positive")
    source_count = 0
    accepted_count = 0
    rejections: Counter[str] = Counter()
    templates: dict[str, GenericCoreTemplate] = {}
    support_units: dict[str, set[str]] = {}
    for row in rows:
        source_count += 1
        result = compile_generic_templates(row, engine=engine, levels=levels)
        if not result.templates:
            rejections[str(result.rejection_reason or "unknown_rejection")] += 1
            continue
        accepted_count += 1
        for compiled in result.templates:
            current = templates.get(compiled.template_id)
            reference = compiled.precedents[0].reference_id
            support_key = reference or compiled.precedents[0].reaction_id
            support_units.setdefault(compiled.template_id, set()).add(support_key)
            if current is None:
                templates[compiled.template_id] = compiled
                continue
            precedents = current.precedents
            if len(precedents) < max_precedents_per_template:
                precedents = (*precedents, compiled.precedents[0])
            templates[compiled.template_id] = replace(
                current,
                observation_support=current.observation_support + 1,
                precedents=tuple(precedents),
            )
    finalized = tuple(
        replace(
            template,
            independent_reference_support=len(support_units[template.template_id]),
        )
        for template in sorted(templates.values(), key=lambda item: item.template_id)
    )
    return GenericTemplateLibrary(
        templates=finalized,
        source_row_count=source_count,
        accepted_observation_count=accepted_count,
        rejection_counts=dict(sorted(rejections.items())),
        definition={**GENERIC_LIBRARY_DEFINITION, "compiler_engine": engine},
    )


def save_generic_library(
    library: GenericTemplateLibrary,
    destination: str | Path,
) -> None:
    """Write canonical gzip JSON."""

    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        library.to_dict(),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as handle:
            handle.write(payload)


def load_generic_library(source: str | Path) -> GenericTemplateLibrary:
    """Load a generic gzip JSON library."""

    with gzip.open(Path(source), "rt", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("generic library must contain an object")
    return GenericTemplateLibrary.from_dict(value)


__all__ = [
    "GENERIC_LIBRARY_DEFINITION",
    "build_generic_library",
    "load_generic_library",
    "save_generic_library",
]
