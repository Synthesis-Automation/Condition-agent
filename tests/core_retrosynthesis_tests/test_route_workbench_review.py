"""SVG review regressions for provider-backed workbench reports."""

from __future__ import annotations

import json

from core_retrosynthesis.route_workbench_review import (
    render_provider_route_workbench_review_html,
    write_provider_route_workbench_review_html,
)


def _report() -> dict:
    step = {
        "product_smiles": "CCN",
        "precursor_smiles": ["CCBr", "N"],
    }
    evidence = {
        "provider_id": "provider:test",
        "provider_rank": 1,
        "operator_id": "operator:test",
        "strategic_class": "scaffold_split",
        "reaction_compatibility_disposition": "pass",
        "selectivity_warning_count": 0,
        "precedent_evidence": {"matches": []},
    }
    return {
        "summary": {"completed_count": 1},
        "cases": [
            {
                "execution_status": "completed",
                "target_name": "Ethylamine <test>",
                "target_smiles": "CCN",
                "challenge_focus": "escaping & SVG",
                "report": {
                    "route_kind": "partial",
                    "routes": [
                        {
                            "route": {
                                "steps": [step],
                                "leaves": [
                                    {
                                        "smiles": "CCBr",
                                        "terminal": False,
                                        "unresolved_reason": "maximum_depth",
                                    },
                                    {
                                        "smiles": "N",
                                        "terminal": True,
                                        "unresolved_reason": None,
                                    },
                                ],
                            },
                            "verification": {
                                "status": "failed",
                                "gates": [
                                    {
                                        "gate": "terminal_completion",
                                        "status": "fail",
                                    }
                                ],
                            },
                            "step_evidence": [evidence],
                            "issues": [
                                {
                                    "issue_id": "issue:1",
                                    "kind": "unresolved_leaf",
                                    "severity": "strong",
                                    "message": "Leaf remains unresolved.",
                                }
                            ],
                            "weakest_issue_id": "issue:1",
                        }
                    ],
                },
            }
        ],
    }


def test_review_renders_inline_target_route_and_leaf_svgs() -> None:
    document = render_provider_route_workbench_review_html(_report())

    assert document.startswith("<!doctype html>")
    assert document.count("<svg") >= 4
    assert "Ethylamine &lt;test&gt;" in document
    assert "escaping &amp; SVG" in document
    assert "Weakest issue" in document
    assert "CCBr.N&gt;&gt;CCN" not in document


def test_review_writer_creates_self_contained_html(tmp_path) -> None:
    source = tmp_path / "report.json"
    output = tmp_path / "review.html"
    source.write_text(json.dumps(_report()), encoding="utf-8")

    result = write_provider_route_workbench_review_html(source, output)

    assert result["case_count"] == 1
    assert result["html_bytes"] == output.stat().st_size
    assert "<svg" in output.read_text(encoding="utf-8")
