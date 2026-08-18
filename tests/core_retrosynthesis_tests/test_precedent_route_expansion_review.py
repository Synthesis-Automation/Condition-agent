"""Self-contained precedent-route expansion review regressions."""

from __future__ import annotations

import json

import pytest

from core_retrosynthesis.cli import main
from core_retrosynthesis.precedent_route_expansion_review import (
    load_precedent_route_expansion_report,
    render_precedent_route_expansion_html,
    write_precedent_route_expansion_html,
)


def _report() -> dict[str, object]:
    pathway = {
        "pathway_id": "path-1",
        "first_abstraction_level": "L0",
        "second_abstraction_level": "L0",
        "first_step_starting_materials": "CCO",
        "intermediate_smiles": "CC=O",
        "second_step_partner_smiles": "N",
        "second_step_starting_materials": "CC=O.N",
        "product_smiles": "CC=N",
        "warnings": [],
    }
    input_evidence = [
        {
            "input_smiles": "CCO",
            "label": "observed alcohol",
            "minimum_level": "R0",
            "role": "first_step",
            "source_kind": "observed",
            "stock_evidence_complete": True,
            "stock_components": [
                {
                    "smiles": "CCO",
                    "canonical_smiles": "CCO",
                    "source_records": [
                        {"supplier": "Example", "snapshot_date": "2026-08-01"}
                    ],
                }
            ],
        }
    ]
    levels = []
    for level in ("R0", "R1", "R2"):
        levels.append(
            {
                "level": level,
                "product_count": 1,
                "product_smiles": ["CC=N"],
                "pathways": [pathway],
                "input_evidence": input_evidence,
            }
        )
    return {
        "artifact_type": "two_step_precedent_route_expansion_poc",
        "algorithm_version": "test@1.0",
        "definition_id": "review-test@1.0",
        "route_count": 1,
        "level_product_counts": {"R0": 1, "R1": 1, "R2": 1},
        "source_validation": {"valid": True, "source_sha256": "abc123"},
        "stock_index_metadata": {
            "source_summaries": [
                {"supplier": "Example", "snapshot_date": "2026-08-01"}
            ]
        },
        "routes": [
            {
                "route_id": "observed_oxidation_condensation",
                "source_route_id": "US_TEST_1",
                "patent_id": "US_TEST",
                "first_source_reaction_id": "rxn-1",
                "second_source_reaction_id": "rxn-2",
                "first_reaction_smiles": "CCO>>CC=O",
                "second_reaction_smiles": "CC=O.N>>CC=N",
                "levels": levels,
            }
        ],
    }


def test_render_review_contains_products_routes_and_review_controls() -> None:
    document = render_precedent_route_expansion_html(_report())

    assert "Two-step precedent-route product review" in document
    assert "Observed two-step precedent" in document
    assert "US_TEST_1" in document
    assert "CC=N" in document
    assert "Pathway 1" in document
    assert "R0 exact" in document
    assert "Export review JSON" in document
    assert "localStorage" in document
    assert 'class="review-step1"' in document
    assert 'class="review-step2"' in document
    assert 'class="render-error"' not in document


def test_write_review_is_deterministic_and_cli_reports_summary(
    tmp_path, capsys
) -> None:
    source = tmp_path / "expansion.json"
    source.write_text(json.dumps(_report()), encoding="utf-8")
    first = tmp_path / "first.html"
    second = tmp_path / "second.html"

    direct_summary = write_precedent_route_expansion_html(source, first)
    assert direct_summary["route_count"] == 1
    assert direct_summary["product_count"] == 1
    first_document = first.read_text(encoding="utf-8")
    write_precedent_route_expansion_html(source, first)
    assert first.read_text(encoding="utf-8") == first_document

    assert main(
        [
            "render-precedent-route-expansion",
            str(source),
            str(second),
            "--title",
            "Review title",
        ]
    ) == 0
    cli_summary = json.loads(capsys.readouterr().out)
    assert cli_summary["product_count"] == 1
    assert "Review title" in second.read_text(encoding="utf-8")


def test_loader_rejects_unrelated_json(tmp_path) -> None:
    source = tmp_path / "unrelated.json"
    source.write_text('{"artifact_type": "other", "routes": [{}]}', encoding="utf-8")

    with pytest.raises(ValueError, match="unexpected.*artifact type"):
        load_precedent_route_expansion_report(source)
