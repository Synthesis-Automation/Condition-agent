from chem_assistant.rule_builder_session import RuleBuilderSession
from chemtools.rule.builder import RuleBuilder


def test_rule_builder_session_wizard_collects_minimum():
    builder = RuleBuilder.new("Suzuki_Miyaura")
    answers = iter(
        [
            # Metadata
            "test_rule",
            "Test Rule",
            "v1",
            "2025-01-01",
            "draft",
            "suzuki,demo",
            # Reference reactions
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
            "",
            # Scope & notes
            "",
            "",
            "",
            "",
            # Applies-if
            "sp2_halide_present,sp2_boron_present",
            "",
            # Default rule
            "DEF",
            "Baseline",
            "pd_source=PdCl2(dtbpf)",
            "ligand=dtbpf",
            "",
            # Base rules
            "y",
            "BR1",
            "Primary",
            "desc",
            "sp2_halide_present",
            "",
            "",
            "pd_source=PdCl2(dtbpf)",
            "ligand=dtbpf",
            "",
            "n",
            # Modifiers
            "y",
            "MOD1",
            "symptom:low_yield",
            "Increase temperature",
            "",
            "n",
        ]
    )

    def fake_input(prompt: str) -> str:
        return next(answers)

    session = RuleBuilderSession(
        builder,
        prompt_fn=fake_input,
        output_fn=lambda *_args, **_kwargs: None,
    )
    session.run_wizard()
    # Should not raise
    builder.validate(strict=True)
    data = builder.to_dict()
    assert data["metadata"]["id"] == "test_rule"
    assert data["base_rules"], "Expected at least one base rule"
    assert data["modifiers"], "Expected at least one modifier"
