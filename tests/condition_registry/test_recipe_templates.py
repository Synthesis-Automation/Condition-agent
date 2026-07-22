import copy
import json

from condition_registry import (
    get_recipe_template,
    load_recipe_template_set,
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
    assert len(template_set.templates) == 5
    assert all(template.identity_complete for template in template_set.templates)
    assert all(template.status == "draft" for template in template_set.templates)
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
        "BH_CN_",
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
