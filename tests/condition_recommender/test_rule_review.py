import csv
import sys

import pytest

from condition_recommender.rule_review import (
    RULE_REVIEW_COLUMNS,
    build_rule_review_rows,
    export_rule_review_csv,
)
from condition_recommender.rule_review_cli import main as rule_review_cli_main
from condition_recommender.rules import load_condition_rule_set
from condition_registry import load_recipe_template_set


def test_rule_review_has_one_row_per_explicit_referenced_variant() -> None:
    rule_set = load_condition_rule_set()
    template_set = load_recipe_template_set()
    templates = {
        template.template_id: template for template in template_set.templates
    }
    expected = sum(
        max(1, len(templates[recommendation.recipe_template_id].variants))
        for rule in rule_set.rules
        for recommendation in rule.recommendations
    )

    rows = build_rule_review_rows()

    assert len(rows) == expected == 15
    assert len(
        {
            (
                row["rule_id"],
                row["recipe_variant_id"],
            )
            for row in rows
        }
    ) == len(rows)
    assert all(tuple(row) == RULE_REVIEW_COLUMNS for row in rows)


def test_rule_review_resolves_conditions_for_chemist_review() -> None:
    rows = build_rule_review_rows()
    active = next(row for row in rows if row["status"] == "active")

    assert active["rule_id"] == "pd_sp2_cn.amide_nh.v1"
    assert active["recipe_variant_id"] == (
        "tbu_brettphos_pd_k3po4_tbuoh.v1"
    )
    assert active["catalyst"] == (
        "t-BuBrettPhos Palladacycle Gen. 3 (1 mol%)"
    )
    assert "Tripotassium phosphate" in active["base"]
    assert active["process_conditions"] == "110 C | 1.5 h | 0.5 M | Ar"
    assert active["stoichiometry"] == (
        "electrophile 1 equiv | nucleophile 1.05-1.2 equiv"
    )
    assert "DOI: 10.1021/ol401208t" in active["source"]
    assert "rule_definition_id" not in active
    assert "recipe_id" not in active


def test_rule_review_summarizes_constraints_and_has_blank_review_fields() -> None:
    rows = build_rule_review_rows()
    hindered = next(
        row
        for row in rows
        if row["rule_id"] == "pd_sp2_cn.free_amine.hindered_ar_cl.v1"
    )

    assert hindered["rule_kind"] == "structural override"
    assert "electrophile: leaving_group" in hindered["match_summary"]
    assert "Cl" in hindered["match_summary"]
    assert "ortho occupancy >= 2" in hindered["match_summary"]
    assert "high" in hindered["match_summary"]
    assert "primary_amine, secondary_amine" in hindered["match_summary"]
    assert hindered["review_decision"] == ""
    assert hindered["review_notes"] == ""


def test_export_rule_review_csv_is_excel_friendly_and_source_consistent(
    tmp_path,
) -> None:
    output = tmp_path / "review" / "pd_sp2_cn.csv"

    built_rows = export_rule_review_csv(output)

    assert output.read_bytes().startswith(b"\xef\xbb\xbf")
    with output.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        exported_rows = list(reader)
    assert tuple(reader.fieldnames or ()) == RULE_REVIEW_COLUMNS
    assert len(RULE_REVIEW_COLUMNS) == 20
    assert exported_rows == list(built_rows)


def test_rule_review_cli_exports_rows(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    output = tmp_path / "rulebook.csv"
    monkeypatch.setattr(
        sys,
        "argv",
        ["rule_review_cli", "export", str(output)],
    )

    rule_review_cli_main()

    assert output.exists()
    assert "Wrote 15 rule-variant rows with 20 columns" in capsys.readouterr().out
