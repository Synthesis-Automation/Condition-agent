"""Render representative sample reactions with the display-minimization POC.

Source labels select review examples only. They are not passed to reaction
featurization, atom mapping, graph minimization, or rendering.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
import re
import sys
from typing import Any, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import (  # noqa: E402
    RxnMapperProvider,
    analyze_reaction_with_external_mapping,
    build_reaction_display_projection,
    featurize_reaction,
    reaction_render_context_from_analysis,
)
from visualization import build_reaction_display_graphic  # noqa: E402


DEFAULT_SOURCE = (
    PROJECT_ROOT / "raw_dataset" / "examples" / "sample_reactions.csv"
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "reaction_display_poc"
DEFAULT_REPORT = DEFAULT_OUTPUT_DIR / "reaction_display_poc.json"

SAMPLE_SELECTORS = (
    ("Suzuki", "Electron-rich ArBr"),
    ("Hydrogenation", "Styrene ethylbenzene"),
    ("Click", "Ether-substituted triazole"),
    ("Amide: Benzoic acid + aniline", ""),
    ("Buchwald-Hartwig amination", "Ph-Br + diethylamine"),
)


def _read_rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return tuple(dict(row) for row in csv.DictReader(handle))


def _selected_rows(
    rows: Sequence[Mapping[str, str]],
) -> tuple[dict[str, str], ...]:
    selected = []
    for reaction_type, example in SAMPLE_SELECTORS:
        matches = [
            row
            for row in rows
            if row.get("reaction_type") == reaction_type
            and row.get("example") == example
        ]
        if len(matches) != 1:
            raise ValueError(
                "expected exactly one sample row for "
                f"{reaction_type!r}/{example!r}; found {len(matches)}"
            )
        selected.append(dict(matches[0]))
    return tuple(selected)


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")


def _analysis_with_core(reaction_smiles: str) -> tuple[Any, dict[str, Any]]:
    baseline = featurize_reaction(reaction_smiles)
    if baseline.reaction_core is not None:
        return baseline, {
            "attempted": False,
            "status": "not_needed",
            "provider": None,
            "confidence": None,
            "warnings": [],
        }
    assessment = analyze_reaction_with_external_mapping(
        reaction_smiles,
        RxnMapperProvider(),
        base_analysis=baseline,
        force_resolved_shadow=True,
    )
    mapping = assessment.mapping_result
    if assessment.analysis.reaction_core is None:
        raise RuntimeError("external mapping did not produce a reaction core")
    return assessment.analysis, {
        "attempted": True,
        "status": assessment.status,
        "provider": assessment.provider_metadata.provider_id,
        "confidence": (
            float(mapping.mapper_confidence) if mapping is not None else None
        ),
        "warnings": list(assessment.warnings),
    }


def build_poc(
    *,
    source_path: Path,
    output_dir: Path,
    report_path: Path,
) -> dict[str, Any]:
    """Build SVGs and a JSON review report for the selected sample rows."""
    rows = _selected_rows(_read_rows(source_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    examples = []
    for row in rows:
        reaction_type = str(row["reaction_type"])
        example = str(row["example"])
        reaction_smiles = str(row["rxn_smiles"])
        analysis, mapping = _analysis_with_core(reaction_smiles)
        projection = build_reaction_display_projection(
            reaction_render_context_from_analysis(analysis)
        )
        graphic = build_reaction_display_graphic(
            projection,
            size=(1200, 300),
        )
        png = build_reaction_display_graphic(
            projection,
            size=(1200, 300),
            image_format="png",
        )
        stem = f"{_slug(reaction_type)}_{_slug(example)}"
        svg_path = output_dir / f"{stem}.svg"
        png_path = output_dir / f"{stem}.png"
        svg_path.write_bytes(graphic.image_bytes)
        png_path.write_bytes(png.image_bytes)
        examples.append(
            {
                "source_selection": {
                    "reaction_type": reaction_type,
                    "example": example,
                    "labels_used_for_chemistry": False,
                },
                "input_reaction_smiles": reaction_smiles,
                "minimum_reaction_smiles": (
                    projection.minimum_reaction_smiles
                ),
                "render_reaction_smiles": projection.render_reaction_smiles,
                "projection": {
                    "schema_version": projection.schema_version,
                    "definition_id": projection.definition_id,
                    "reaction_scope": projection.reaction_scope,
                    "formed_ring_sizes": list(projection.formed_ring_sizes),
                    "annotation": projection.annotation,
                    "evidence_status": projection.evidence_status,
                    "confidence": projection.confidence,
                    "warnings": list(projection.warnings),
                    "reactants": [
                        asdict(component) for component in projection.reactants
                    ],
                    "products": [
                        asdict(component) for component in projection.products
                    ],
                },
                "mapping": mapping,
                "svg_path": svg_path.relative_to(PROJECT_ROOT).as_posix(),
                "png_path": png_path.relative_to(PROJECT_ROOT).as_posix(),
            }
        )
    report = {
        "schema_version": "1.0",
        "artifact_type": "reaction_display_minimization_poc",
        "source_csv": source_path.relative_to(PROJECT_ROOT).as_posix(),
        "source_labels_used_for_chemistry": False,
        "placeholder_smiles_meaning": {"*": "R"},
        "example_count": len(examples),
        "examples": examples,
        "limitations": [
            "Display projections are not reaction identities.",
            "Unchanged substituents removed from aromatic carbon remain "
            "available in the complete reaction analysis.",
            "An external atom-mapping proposal is review evidence until "
            "independently validated.",
        ],
    }
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args(argv)
    report = build_poc(
        source_path=args.source,
        output_dir=args.output_dir,
        report_path=args.report,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
