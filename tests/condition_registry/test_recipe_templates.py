import copy
import json

import pytest

from condition_registry import (
    get_recipe_template,
    load_recipe_template_set,
    materialize_recipe_variant,
    resolve_substance_id,
)
from condition_registry.template_loader import (
    RECIPE_TEMPLATES_PATH,
    validate_recipe_template_payload,
)


def _payload() -> dict:
    return json.loads(RECIPE_TEMPLATES_PATH.read_text(encoding="utf-8"))


def test_clean_recipe_templates_are_typed_and_registry_backed() -> None:
    template_set = load_recipe_template_set()

    assert template_set.definition_id == "condition_recipe_templates.v1"
    assert template_set.schema_version == "1.2"
    assert len(template_set.templates) == 11
    assert all(template.identity_complete for template in template_set.templates)
    statuses = {template.template_id: template.status for template in template_set.templates}
    assert statuses["pd_sp2_cn.amide_nh.v1"] == "active"
    assert sum(status == "active" for status in statuses.values()) == 1
    for template in template_set.templates:
        for slot in template.slots:
            for option in slot.alternatives:
                result = resolve_substance_id(option.substance_id)
                assert result.status == "resolved"
                assert result.substance is not None
                assert slot.role_id in {
                    assignment.role_id for assignment in result.substance.roles
                }


def test_recipe_templates_do_not_embed_legacy_rule_fields() -> None:
    serialized = json.dumps(_payload(), sort_keys=True)

    for legacy_token in (
        "reactant_1",
        "reactant_2",
        "rule_tier",
        "z_score",
        "datasets/rules",
    ):
        assert legacy_token not in serialized


def test_recipe_alternatives_are_not_implicitly_materialized() -> None:
    template = get_recipe_template("pd_sp2_cn.primary_alkyl_amine.v1")

    assert template is not None
    base = next(slot for slot in template.slots if slot.slot_id == "base")
    assert base.selection_policy == "present_alternatives"
    assert len(base.alternatives) == 4


def test_template_validation_rejects_unknown_and_role_mismatched_identities() -> None:
    payload = _payload()
    unknown = copy.deepcopy(payload)
    unknown["templates"][0]["slots"][0]["alternatives"][0][
        "substance_id"
    ] = "cas:0000-00-0"
    assert any(
        "unknown_substance_id" in error
        for error in validate_recipe_template_payload(unknown)
    )

    mismatched = copy.deepcopy(payload)
    mismatched["templates"][0]["slots"][0]["alternatives"][0][
        "substance_id"
    ] = "cas:75-05-8"
    assert any(
        "substance_role_mismatch" in error
        for error in validate_recipe_template_payload(mismatched)
    )


def test_template_validation_rejects_legacy_or_unknown_fields() -> None:
    payload = _payload()
    payload["templates"][0]["yield"] = 95

    assert "templates[0]:unknown_key:yield" in validate_recipe_template_payload(
        payload
    )


def test_explicit_variant_materializes_one_canonical_recipe() -> None:
    template = get_recipe_template("pd_sp2_cn.primary_alkyl_amine.v1")
    assert template is not None

    recipe = materialize_recipe_variant(
        template,
        "tbu_brettphos_k3po4_mecn.v1",
        transformation_class="sp2_c_n_substitution",
        include_draft=True,
    )

    assert recipe.recipe_id.startswith("RCR1:")
    assert len(recipe.components) == 3
    assert recipe.catalysts[0].substance_id == "cas:1536473-72-9"
    assert recipe.bases[0].substance_id == "cas:7778-53-2"
    assert recipe.solvents[0].substance_id == "cas:75-05-8"
    assert all(
        component.provenance["recipe_variant_id"]
        == "tbu_brettphos_k3po4_mecn.v1"
        for component in recipe.components
    )
    assert len(template.variants) == 2
    assert len(next(slot for slot in template.slots if slot.slot_id == "base").alternatives) == 4

    repeated = materialize_recipe_variant(
        template,
        "tbu_brettphos_k3po4_mecn.v1",
        transformation_class="sp2_c_n_substitution",
        include_draft=True,
    )
    assert recipe.recipe_id == repeated.recipe_id
    assert recipe.catalysts[0].amount == 2.0
    assert recipe.catalysts[0].amount_unit == "mol_percent"
    assert recipe.bases[0].amount == 2.0
    assert recipe.bases[0].amount_unit == "equivalent"
    assert recipe.time_h == 16.0
    assert recipe.concentration_m == 0.2


def test_every_draft_screening_variant_is_explicitly_materializable() -> None:
    template_set = load_recipe_template_set()

    for template in template_set.templates:
        if template.status != "draft":
            continue
        assert template.temperature_c is not None
        assert template.time_h is not None
        assert template.concentration_m is not None
        assert template.atmosphere
        assert {amount.role for amount in template.partner_amounts} == {
            "electrophile",
            "nucleophile",
        }
        for variant in template.variants:
            recipe = materialize_recipe_variant(
                template,
                variant.variant_id,
                transformation_class="sp2_c_n_substitution",
                include_draft=True,
            )
            assert recipe.catalysts
            assert recipe.bases
            assert recipe.solvents
            assert all(
                component.amount is not None
                for component in (*recipe.catalysts, *recipe.bases)
            )


def test_draft_variant_materialization_requires_explicit_opt_in() -> None:
    template = get_recipe_template("pd_sp2_cn.primary_alkyl_amine.v1")
    assert template is not None

    with pytest.raises(ValueError, match="explicit opt-in"):
        materialize_recipe_variant(
            template,
            "tbu_brettphos_k3po4_mecn.v1",
            transformation_class="sp2_c_n_substitution",
        )


def test_template_validation_rejects_implicit_or_invalid_variant_choices() -> None:
    payload = _payload()
    variant = payload["templates"][0]["variants"][0]
    base = next(
        selection
        for selection in variant["selections"]
        if selection["slot_id"] == "base"
    )
    base["substance_id"] = "cas:123-91-1"

    assert any(
        "selection_not_in_slot_alternatives:base:cas:123-91-1" in error
        for error in validate_recipe_template_payload(payload)
    )


def test_active_template_is_complete_and_materializes_quantities() -> None:
    template = get_recipe_template("pd_sp2_cn.amide_nh.v1")
    assert template is not None

    recipe = materialize_recipe_variant(
        template,
        "tbu_brettphos_pd_k3po4_tbuoh.v1",
        transformation_class="sp2_c_n_substitution",
    )

    assert recipe.temperature_c == 110.0
    assert recipe.time_h == 1.5
    assert recipe.concentration_m == 0.5
    assert recipe.atmosphere == "Ar"
    assert recipe.catalysts[0].amount == 1.0
    assert recipe.catalysts[0].amount_unit == "mol_percent"
    assert recipe.bases[0].amount == 1.4
    assert recipe.bases[0].amount_unit == "equivalent"
    assert recipe.solvents[0].substance_id == "cas:75-65-0"
    assert tuple(amount.role for amount in template.partner_amounts) == (
        "electrophile",
        "nucleophile",
    )


def test_activation_validation_rejects_incomplete_production_protocol() -> None:
    payload = _payload()
    amide = next(
        template
        for template in payload["templates"]
        if template["template_id"] == "pd_sp2_cn.amide_nh.v1"
    )
    del amide["time_h"]
    amide["variants"][0]["selections"][0].pop("amount")
    amide["variants"][0]["selections"][0].pop("amount_unit")
    amide["provenance"]["sources"] = []

    errors = validate_recipe_template_payload(payload)

    assert any("active_template_missing:time_h" in error for error in errors)
    assert any("active_variant_missing_quantity:metal_catalyst" in error for error in errors)
    assert any("active_template_missing:provenance_sources" in error for error in errors)
