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
                row["recipe_template_id"],
                row["recipe_variant_id"],
            )
            for row in rows
        }
    ) == len(rows)
    assert all(tuple(row) == RULE_REVIEW_COLUMNS for row in rows)


def test_rule_review_resolves_conditions_and_preserves_definition_ids() -> None:
    rows = build_rule_review_rows()
    active = next(row for row in rows if row["rule_status"] == "active")

    assert active["rule_id"] == "pd_sp2_cn.amide_nh.v1"
    assert active["recipe_variant_id"] == (
        "tbu_brettphos_pd_k3po4_tbuoh.v1"
    )
    assert (
        "t-BuBrettPhos Palladacycle Gen. 3 "
        "[cas:1536473-72-9] (1 mol_percent)"
        == active["catalyst"]
    )
    assert "Tripotassium phosphate" in active["base"]
    assert active["temperature_c"] == "110.0"
    assert active["nucleophile_equiv_min"] == "1.05"
    assert active["doi"] == "10.1021/ol401208t"
    assert active["recipe_id"].startswith("RCR1:")
    assert active["rule_definition_id"] == "pd_sp2_cn_condition_rules.v1"
    assert active["template_definition_id"] == (
        "condition_recipe_templates.v1"
    )


def test_rule_review_exposes_structural_constraints_and_blank_review_fields() -> None:
    rows = build_rule_review_rows()
    hindered = next(
        row
        for row in rows
        if row["rule_id"] == "pd_sp2_cn.free_amine.hindered_ar_cl.v1"
    )

    assert hindered["rule_kind"] == "structural_override"
    assert hindered["electrophile_handle_tokens"] == "Cl"
    assert hindered["electrophile_steric_classes"] == "ortho_hindered"
    assert hindered["electrophile_ortho_substituent_count_min"] == "2"
    assert hindered["nucleophile_families"] == (
        "primary_amine | secondary_amine"
    )
    assert hindered["legacy_source_rule_ids"] == (
        "BH_CN_01 | BH_CN_02 | BH_CN_05"
    )
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
    assert "Wrote 15 rule-variant rows" in capsys.readouterr().out
