"""Render the complex click example as a minimized graphical reaction."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
from typing import Any, Sequence


PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from reactive_taxonomy import (  # noqa: E402
    RxnMapperProvider,
    analyze_reaction_with_external_mapping,
    featurize_reaction,
)
from visualization import build_reaction_core_graphic  # noqa: E402


CLICK_REACTION = (
    "C#CCN(c1nc2ccc(F)cc2s1)c1nc2ccc(OC)cc2s1."
    "[N-]=[N+]=NCC(=O)Nc1ccc(SC(F)(F)F)cc1"
    ">>COc1ccc2nc(N(Cc3cn(CC(=O)Nc4ccc(SC(F)(F)F)cc4)nn3)"
    "c3nc4ccc(F)cc4s3)sc2c1"
)
DEFAULT_SVG = PROJECT_ROOT / "results" / "click_minimized_reaction_poc.svg"
DEFAULT_PNG = PROJECT_ROOT / "results" / "click_minimized_reaction_poc.png"
DEFAULT_REPORT = (
    PROJECT_ROOT / "results" / "click_minimized_reaction_poc.json"
)


def _placeholder_payload(value: Any) -> dict[str, Any]:
    return {
        "label": value.label,
        "remote_class": value.remote_class,
        "fragment_smiles": value.fragment_smiles,
        "functional_group_ids": list(value.functional_group_ids),
    }


def build_poc(
    *,
    svg_path: Path,
    png_path: Path,
    report_path: Path,
) -> dict[str, Any]:
    """Map, minimize, render, and persist the click-reaction POC."""
    base = featurize_reaction(CLICK_REACTION)
    assessment = analyze_reaction_with_external_mapping(
        CLICK_REACTION,
        RxnMapperProvider(),
        base_analysis=base,
        force_resolved_shadow=True,
    )
    analysis = assessment.analysis
    core = analysis.reaction_core
    if core is None:
        raise RuntimeError("RXNMapper did not produce a drawable reaction core")
    svg = build_reaction_core_graphic(
        analysis,
        size=(1200, 260),
        image_format="svg",
    )
    png = build_reaction_core_graphic(
        analysis,
        size=(1200, 260),
        image_format="png",
    )
    for path in (svg_path, png_path, report_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.write_bytes(svg.image_bytes)
    png_path.write_bytes(png.image_bytes)
    mapping = assessment.mapping_result
    report = {
        "schema_version": "1.0",
        "artifact_type": "reaction_core_graphic_poc",
        "expected_reaction_type": "click",
        "source_label_used_for_analysis": False,
        "input_reaction_smiles": CLICK_REACTION,
        "base_analysis": {
            "evidence_quality": base.evidence_quality,
            "has_signature": base.reaction_signature is not None,
            "has_reaction_core": base.reaction_core is not None,
            "named_family": base.named_family,
        },
        "external_mapping": {
            "status": assessment.status,
            "provider": assessment.provider_metadata.provider_id,
            "confidence": (
                mapping.mapper_confidence if mapping is not None else None
            ),
            "mapped_reaction_smiles": (
                mapping.mapped_reaction_smiles if mapping is not None else None
            ),
            "warnings": list(assessment.warnings),
        },
        "reaction_core": {
            "core_id": core.core_id,
            "shape_core_key": core.shape_core_key,
            "center_transition_key": core.center_transition_key,
            "generic_label": core.generic_label,
            "active_atom_count": core.active_atom_count,
            "event_count": core.event_count,
            "evidence_status": core.evidence_status,
            "confidence": core.confidence,
            "participant_tokens": list(core.participant_tokens),
            "warnings": list(core.warnings),
            "atom_transitions": [
                asdict(transition) for transition in core.atom_transitions
            ],
        },
        "graphic": {
            "definition_id": svg.definition_id,
            "schema_version": svg.schema_version,
            "placeholders": [
                _placeholder_payload(value) for value in svg.placeholders
            ],
            "svg_path": str(svg_path),
            "png_path": str(png_path),
        },
        "feasibility": {
            "drawable": True,
            "active_graph": (
                "two alkyne carbons plus three azide nitrogens form a "
                "five-membered triazole"
            ),
            "abstraction": (
                "two unchanged retained substituents are rendered as R1 and R2"
            ),
            "limitation": (
                "low-confidence external mapping makes this a review graphic, "
                "not verified reaction identity"
            ),
        },
    }
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--svg", type=Path, default=DEFAULT_SVG)
    parser.add_argument("--png", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args(argv)
    report = build_poc(
        svg_path=args.svg,
        png_path=args.png,
        report_path=args.report,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
