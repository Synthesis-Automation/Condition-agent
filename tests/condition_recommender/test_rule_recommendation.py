import json
import sys

import pytest

from condition_recommender import recommend_rule_conditions
from condition_recommender.compatibility import CompatibilityAssessment
from condition_recommender.rules import api as rule_api
from reactive_taxonomy import featurize_reaction
from condition_recommender.rule_recommend_cli import main as rule_cli_main


@pytest.mark.parametrize(
    ("reaction", "expected_rule"),
    (
        (
            "Brc1ccccc1.CN>>CNc1ccccc1",
            "pd_sp2_cn.primary_alkyl_amine.v1",
        ),
        (
            "Brc1ccccc1.CNC>>CN(C)c1ccccc1",
            "pd_sp2_cn.secondary_alkyl_amine.v1",
        ),
        (
            "Brc1ccccc1.c1cc[nH]c1>>c1ccc(-n2cccc2)cc1",
            "pd_sp2_cn.aromatic_nh.v1",
        ),
        (
            "Clc1ccccc1.CC(N)=O>>CC(=O)Nc1ccccc1",
            "pd_sp2_cn.amide_nh.v1",
        ),
        (
            "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
            "pd_sp2_cn.primary_aryl_amine.v1",
        ),
        (
            "Brc1ccccc1.CNc1ccccc1>>CN(c1ccccc1)c1ccccc1",
            "pd_sp2_cn.secondary_aryl_amine.v1",
        ),
    ),
)
def test_specific_structural_rules_match_without_named_family(
    reaction: str, expected_rule: str
) -> None:
    analysis = featurize_reaction(reaction)
    result = recommend_rule_conditions(reaction, include_draft=True)

    assert analysis.named_family is None
    assert result.valid
    assert result.reaction_signature_schema_version == "3.0"
    assert result.schema_version == "2.0"
    assert result.rule_definition_schema_version == "2.0"
    assert result.recipe_template_schema_version == "1.2"
    assert result.transformation_class == "sp2_c_n_substitution"
    assert dict(result.taxonomy_definition_versions)[
        "reaction_grammars.v2.json"
    ].startswith("2.0@sha256:")
    assert result.selected_tiers == (
        ("pd_sp2_cn.condition_regime", "specific"),
    )
    assert tuple(item.rule_id for item in result.recommendations) == (
        expected_rule,
    )
    if expected_rule == "pd_sp2_cn.amide_nh.v1":
        assert "DRAFT_RULES_INCLUDED" not in result.warnings
    else:
        assert "DRAFT_RULES_INCLUDED" in result.warnings


def test_specific_rule_suppresses_matching_fallback() -> None:
    reaction = "Brc1ccccc1.CN>>CNc1ccccc1"
    result = recommend_rule_conditions(reaction, include_draft=True)
    traces = {trace.rule_id: trace for trace in result.match_traces}

    assert traces["pd_sp2_cn.primary_alkyl_amine.v1"].matched
    assert traces["pd_sp2_cn.general_nh_fallback.v1"].matched
    assert tuple(item.rule_id for item in result.recommendations) == (
        "pd_sp2_cn.primary_alkyl_amine.v1",
    )
    assert result.recommendations[0].rule_kind == "default"


@pytest.mark.parametrize(
    ("reaction", "expected_rule", "temperature_c", "time_h"),
    (
        (
            "Clc1ccccc1.CN>>CNc1ccccc1",
            "pd_sp2_cn.free_amine.ar_cl.v1",
            100.0,
            16.0,
        ),
        (
            "Cc1cccc(C)c1Br.CN>>Cc1cccc(C)c1NC",
            "pd_sp2_cn.free_amine.hindered_ar_br.v1",
            110.0,
            18.0,
        ),
        (
            "Cc1cccc(C)c1Cl.CN>>Cc1cccc(C)c1NC",
            "pd_sp2_cn.free_amine.hindered_ar_cl.v1",
            110.0,
            24.0,
        ),
        (
            "Brc1ccccc1.CC(C)(C)N>>CC(C)(C)Nc1ccccc1",
            "pd_sp2_cn.free_amine.hindered_nucleophile_ar_br.v1",
            110.0,
            18.0,
        ),
    ),
)
def test_structural_overrides_suppress_class_defaults(
    reaction: str,
    expected_rule: str,
    temperature_c: float,
    time_h: float,
) -> None:
    result = recommend_rule_conditions(reaction, include_draft=True)

    assert tuple(item.rule_id for item in result.recommendations) == (
        expected_rule,
    )
    recommendation = result.recommendations[0]
    assert recommendation.rule_kind == "structural_override"
    assert "STRUCTURAL_OVERRIDE_APPLIED" in result.warnings
    assert "DEFAULT_RULE_APPLIED" not in result.warnings
    variant = recommendation.compatible_variants[0]
    assert variant.resolved_recipe["temperature_c"] == temperature_c
    assert variant.resolved_recipe["time_h"] == time_h


def test_combined_hindered_ar_cl_has_highest_priority() -> None:
    reaction = "Cc1cccc(C)c1Cl.CN>>Cc1cccc(C)c1NC"
    result = recommend_rule_conditions(reaction, include_draft=True)
    traces = {trace.rule_id: trace for trace in result.match_traces}

    assert traces["pd_sp2_cn.free_amine.ar_cl.v1"].matched
    assert traces["pd_sp2_cn.free_amine.hindered_ar_cl.v1"].matched
    assert result.recommendations[0].rule_id == (
        "pd_sp2_cn.free_amine.hindered_ar_cl.v1"
    )


def test_draft_override_blocks_unsafe_production_default_fallback() -> None:
    result = recommend_rule_conditions(
        "Clc1ccccc1.CN>>CNc1ccccc1"
    )

    assert result.valid
    assert result.recommendations == ()
    assert result.selected_tiers == ()
    assert "DRAFT_RULE_MATCHES_EXCLUDED" in result.warnings


def test_unsupported_aryl_fluoride_abstains() -> None:
    result = recommend_rule_conditions(
        "Fc1ccccc1.CN>>CNc1ccccc1",
        include_draft=True,
    )

    assert result.valid
    assert result.recommendations == ()
    assert "NO_STRUCTURAL_RULE_MATCH" in result.warnings


def test_intramolecular_cn_cyclization_is_outside_phase_one_scope() -> None:
    result = recommend_rule_conditions(
        "NCCCc1ccccc1Br>>c1ccc2c(c1)CCCN2",
        include_draft=True,
    )

    assert result.valid
    assert result.transformation_class == "sp2_c_n_substitution"
    assert result.recommendations == ()
    assert "NO_STRUCTURAL_RULE_MATCH" in result.warnings
    assert all(
        "reaction_scope:intramolecular" in trace.failures
        for trace in result.match_traces
    )


def test_aniline_uses_primary_aryl_amine_default() -> None:
    reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"
    result = recommend_rule_conditions(reaction, include_draft=True)
    traces = {trace.rule_id: trace for trace in result.match_traces}

    assert result.selected_tiers == (
        ("pd_sp2_cn.condition_regime", "specific"),
    )
    assert tuple(item.rule_id for item in result.recommendations) == (
        "pd_sp2_cn.primary_aryl_amine.v1",
    )
    assert result.recommendations[0].rule_kind == "default"
    primary = traces["pd_sp2_cn.primary_alkyl_amine.v1"]
    assert not primary.matched
    assert "nucleophile:missing_retained_contexts:Alkyl" in primary.failures
    assert "nucleophile:disallowed_retained_contexts:Ar" in primary.failures


def test_heteroaryl_primary_amine_uses_aryl_amine_default() -> None:
    result = recommend_rule_conditions(
        "Brc1ccccc1.Nc1ccccn1>>c1ccc(Nc2ccccn2)cc1",
        include_draft=True,
    )

    assert result.recommendations[0].rule_id == (
        "pd_sp2_cn.primary_aryl_amine.v1"
    )
    assert any(
        "retained_contexts=HeteroAr" in evidence
        for evidence in result.recommendations[0].match_evidence
    )


def test_draft_rules_are_excluded_by_default() -> None:
    result = recommend_rule_conditions(
        "Brc1ccccc1.CN>>CNc1ccccc1"
    )

    assert result.valid
    assert result.recommendations == ()
    assert result.selected_tiers == ()
    assert "DRAFT_RULE_MATCHES_EXCLUDED" in result.warnings
    assert any(trace.matched and not trace.eligible for trace in result.match_traces)


def test_reviewed_primary_amide_protocol_is_available_in_production() -> None:
    result = recommend_rule_conditions(
        "Clc1ccccc1.CC(N)=O>>CC(=O)Nc1ccccc1"
    )

    assert result.valid
    assert result.warnings == ("DEFAULT_RULE_APPLIED",)
    assert result.selected_tiers == (
        ("pd_sp2_cn.condition_regime", "specific"),
    )
    recommendation = result.recommendations[0]
    assert recommendation.rule_status == "active"
    assert recommendation.recipe_template["status"] == "active"
    assert recommendation.recipe_template["provenance"]["sources"][0][
        "doi"
    ] == "10.1021/ol401208t"
    variant = recommendation.compatible_variants[0]
    assert variant.rank == 1
    assert variant.definition_priority == 100
    assert variant.variant_status == "active"
    assert variant.resolved_recipe["temperature_c"] == 110.0
    assert variant.resolved_recipe["time_h"] == 1.5
    assert variant.resolved_recipe["concentration_m"] == 0.5
    assert variant.resolved_recipe["catalysts"][0]["amount_unit"] == "mol_percent"
    assert variant.resolved_recipe["bases"][0]["amount_unit"] == "equivalent"


def test_primary_amide_protocol_does_not_expand_beyond_aryl_chlorides() -> None:
    result = recommend_rule_conditions(
        "Brc1ccccc1.CC(N)=O>>CC(=O)Nc1ccccc1"
    )
    traces = {trace.rule_id: trace for trace in result.match_traces}

    assert result.recommendations == ()
    assert "NO_STRUCTURAL_RULE_MATCH" in result.warnings
    assert "electrophile:handle_token:Br" in traces[
        "pd_sp2_cn.amide_nh.v1"
    ].failures


def test_primary_amide_protocol_rejects_secondary_amides() -> None:
    result = recommend_rule_conditions(
        "Clc1ccccc1.CC(=O)NC>>CC(=O)N(C)c1ccccc1"
    )
    traces = {trace.rule_id: trace for trace in result.match_traces}

    assert result.recommendations == ()
    assert "nucleophile:h_count_below_minimum:1" in traces[
        "pd_sp2_cn.amide_nh.v1"
    ].failures


def test_rule_matching_is_invariant_to_reactant_order() -> None:
    forward = recommend_rule_conditions(
        "Brc1ccccc1.CN>>CNc1ccccc1", include_draft=True
    )
    reverse = recommend_rule_conditions(
        "CN.Brc1ccccc1>>CNc1ccccc1", include_draft=True
    )

    assert forward.query_signature_id == reverse.query_signature_id
    assert tuple(item.rule_id for item in forward.recommendations) == tuple(
        item.rule_id for item in reverse.recommendations
    )


def test_unrelated_transformation_returns_no_structural_rule_match() -> None:
    result = recommend_rule_conditions(
        "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1",
        include_draft=True,
    )

    assert result.valid
    assert result.transformation_class == "sp2_c_o_substitution"
    assert result.recommendations == ()
    assert "NO_STRUCTURAL_RULE_MATCH" in result.warnings


def test_query_without_verified_signature_is_not_rule_eligible() -> None:
    result = recommend_rule_conditions(
        "Brc1ccccc1.CN(C)C>>c1ccccc1",
        include_draft=True,
    )

    assert not result.valid
    assert result.recommendations == ()
    assert result.error == "QUERY_HAS_NO_USABLE_REACTION_SIGNATURE"


def test_recommendation_contains_auditable_facts_and_template_status() -> None:
    result = recommend_rule_conditions(
        "Brc1ccccc1.CN>>CNc1ccccc1", include_draft=True
    )
    recommendation = result.recommendations[0]

    assert "transformation_class=sp2_c_n_substitution" in (
        recommendation.match_evidence
    )
    assert any(
        "derived_family=primary_amine" in item
        for item in recommendation.match_evidence
    )
    assert recommendation.recipe_template["status"] == "draft"
    assert recommendation.recipe_template["identity_complete"] is True
    assert len(recommendation.compatible_variants) == 2
    assert all(
        variant.recipe_id.startswith("RCR1:")
        for variant in recommendation.compatible_variants
    )
    assert all(
        variant.resolved_recipe["catalysts"]
        and variant.resolved_recipe["bases"]
        and variant.resolved_recipe["solvents"]
        for variant in recommendation.compatible_variants
    )
    assert any("not active for production" in item for item in recommendation.cautions)


def test_template_variants_receive_spectator_compatibility_penalties() -> None:
    result = recommend_rule_conditions(
        "Brc1ccc(C(=O)O)cc1.CN>>CNc1ccc(C(=O)O)cc1",
        include_draft=True,
    )

    variants = result.recommendations[0].compatible_variants
    assert len(variants) == 2
    assert all(variant.compatibility_score == 0.84 for variant in variants)
    assert all(
        variant.penalty_ids
        == ("acidic_group_with_base", "metal_binding_with_catalyst")
        for variant in variants
    )
    assert all(variant.compatibility_evidence for variant in variants)


def test_hard_incompatible_variants_are_not_recommended(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        rule_api,
        "assess_recipe_compatibility",
        lambda signature, recipe: CompatibilityAssessment(
            compatible=False,
            score=0.0,
            hard_conflicts=("test_hard_conflict",),
            evidence=("Test hard conflict",),
        ),
    )

    result = recommend_rule_conditions(
        "Brc1ccccc1.CN>>CNc1ccccc1", include_draft=True
    )

    assert result.recommendations == ()
    assert len(result.excluded_variants) == 2
    assert all(not variant.compatible for variant in result.excluded_variants)
    assert "INCOMPATIBLE_RULE_VARIANTS_EXCLUDED:2" in result.warnings
    assert "DRAFT_RULE_MATCHES_EXCLUDED" not in result.warnings


def test_rule_cli_requires_explicit_draft_opt_in(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rule-recommend",
            "Brc1ccccc1.CN>>CNc1ccccc1",
            "--include-draft",
            "--json",
        ],
    )

    rule_cli_main()

    payload = json.loads(capsys.readouterr().out)
    assert payload["recommendations"][0]["rule_id"] == (
        "pd_sp2_cn.primary_alkyl_amine.v1"
    )
    assert "DRAFT_RULES_INCLUDED" in payload["warnings"]


def test_rule_cli_defaults_to_concise_review_output(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rule-recommend",
            "Clc1ccccc1.CN>>CNc1ccccc1",
            "--include-draft",
        ],
    )

    rule_cli_main()

    output = capsys.readouterr().out
    assert "Transformation: sp2_c_n_substitution" in output
    assert "AdBrettPhos Pd G3, 2 mol%" in output
    assert "Sodium tert-butoxide, 1.5 equiv" in output
    assert "Conditions: 100 C; 16 h; 0.2 M; N2" in output
    assert "Stoichiometry: electrophile 1 equiv; nucleophile 1.1-1.2 equiv" in output
    assert "match_traces" not in output
    assert '"recipe_template"' not in output
