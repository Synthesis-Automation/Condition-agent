"""Template-library construction and canonical JSON persistence."""

from __future__ import annotations

import gzip
import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, TextIO

from ...row_io import iter_rows

from .extraction import ExtractedTemplate, extract_cx_template
from .models import (
    CxTemplate,
    LibraryBuildReport,
    RetrosynthesisLibrary,
    TemplatePrecedent,
)


LIBRARY_DEFINITION = {
    "definition_id": "cx_retrosynthesis_poc.v2",
    "allowed_bonds": ["C-N", "C-O", "C-S"],
    "allowed_core_quality": ["pass", "review"],
    "allowed_core_evidence": ["verified", "inferred"],
    "required_heavy_atom_event_count": 1,
    "attached_hydrogen_loss_events_allowed": True,
    "required_product_count": 1,
    "required_product_completeness": "verified",
    "source_round_trip_required": True,
    "inferred_mapping_must_revalidate_as_pass": True,
    "materialized_mapping_center_must_match_source_core": True,
    "stereo_changes_allowed": False,
    "bond_order_changes_allowed": False,
}


@dataclass
class _Accumulator:
    extracted: ExtractedTemplate
    observation_support: int = 0
    support_units: set[str] = field(default_factory=set)
    precedents: list[TemplatePrecedent] = field(default_factory=list)


def _open_library_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def build_library(
    rows: Iterable[Dict[str, Any]],
    *,
    max_rows: Optional[int] = None,
    max_precedents_per_template: int = 8,
) -> tuple[RetrosynthesisLibrary, LibraryBuildReport]:
    """Extract, round-trip, deduplicate, and aggregate C-X templates."""

    if max_rows is not None and max_rows < 1:
        raise ValueError("max_rows must be positive")
    if max_precedents_per_template < 1:
        raise ValueError("max_precedents_per_template must be positive")
    accumulators: dict[str, _Accumulator] = {}
    rejections: Counter[str] = Counter()
    source_count = 0
    accepted_count = 0
    for row in rows:
        if max_rows is not None and source_count >= max_rows:
            break
        source_count += 1
        result = extract_cx_template(row)
        if result.template is None:
            rejections[str(result.rejection_reason or "unknown_rejection")] += 1
            continue
        accepted_count += 1
        extracted = result.template
        accumulator = accumulators.setdefault(
            extracted.template_id,
            _Accumulator(extracted=extracted),
        )
        accumulator.observation_support += 1
        accumulator.support_units.add(extracted.precedent.support_unit_id)
        if (
            len(accumulator.precedents) < max_precedents_per_template
            and extracted.precedent not in accumulator.precedents
        ):
            accumulator.precedents.append(extracted.precedent)

    templates = []
    for template_id, accumulator in sorted(accumulators.items()):
        extracted = accumulator.extracted
        templates.append(
            CxTemplate(
                template_id=template_id,
                bond_kind=extracted.bond_kind,
                reaction_smarts=extracted.reaction_smarts,
                product_smarts=extracted.product_smarts,
                precursor_smarts=extracted.precursor_smarts,
                intra_only=extracted.intra_only,
                dimer_only=extracted.dimer_only,
                observation_support=accumulator.observation_support,
                independent_reference_support=len(accumulator.support_units),
                precedents=tuple(accumulator.precedents),
            )
        )
    rejection_counts = dict(sorted(rejections.items()))
    library = RetrosynthesisLibrary(
        templates=tuple(templates),
        source_row_count=source_count,
        accepted_observation_count=accepted_count,
        rejection_counts=rejection_counts,
        definition=dict(LIBRARY_DEFINITION),
    )
    report = LibraryBuildReport(
        source_row_count=source_count,
        accepted_observation_count=accepted_count,
        unique_template_count=len(templates),
        rejection_counts=rejection_counts,
    )
    return library, report


def save_library(library: RetrosynthesisLibrary, destination: str | Path) -> None:
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


def load_library(source: str | Path) -> RetrosynthesisLibrary:
    """Load and validate a JSON or JSON.GZ template library."""

    path = Path(source)
    with _open_library_text(path) as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("retrosynthesis library must be a JSON object")
    return RetrosynthesisLibrary.from_dict(value)


__all__ = [
    "LIBRARY_DEFINITION",
    "build_library",
    "iter_rows",
    "load_library",
    "save_library",
]
