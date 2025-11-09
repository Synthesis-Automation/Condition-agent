import pytest

from chemtools.rule.builder import RuleBuilder


def _seed_builder() -> RuleBuilder:
    builder = RuleBuilder.new("Suzuki_Miyaura")
    builder.set_metadata(
        id="test_rule",
        name="Test Rule",
        version="v0",
        created_date="2025-01-01",
        status="draft",
        tags=["suzuki", "demo"],
    )
    builder.add_reference_reactions(
        ["Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"]
    )
    builder.set_applies_if(all_features=["sp2_halide_present", "sp2_boron_present"])
    builder.set_default_rule(
        rule_id="default",
        description="baseline",
        conditions={
            "pd_source": "PdCl2(dtbpf)",
            "ligand": "dtbpf",
            "base": "K3PO4",
            "temperature_C": "100",
            "time_h": "2",
        },
    )
    return builder


def test_builder_happy_path():
    builder = _seed_builder()
    builder.upsert_base_rule(
        "rule_a",
        name="Rule A",
        description="demo",
        reactant_features={"all": ["sp2_halide_present"]},
        conditions={
            "pd_source": "PdCl2(dtbpf)",
            "ligand": "dtbpf",
            "base": "K3PO4",
            "temperature_C": "100",
            "time_h": "2",
        },
    )
    builder.upsert_modifier(
        "mod_a",
        when=["symptom:low_yield"],
        suggestion="Increase temperature",
    )

    issues = builder.validate(strict=False)
    assert issues == []


def test_builder_detects_missing_rules():
    builder = _seed_builder()
    with pytest.raises(ValueError) as exc:
        builder.validate(strict=True)
    assert "base_rules" in str(exc.value)
