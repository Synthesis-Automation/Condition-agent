import json

from chem_assistant.chemtools_wrapper import (
    RuleBuilderAutoInput,
    run_rule_builder_autofill,
)


def test_rule_builder_autofill_success(monkeypatch):
    payload = {
        "notes": "HTE plates highlight dtbpf/K3PO4 dominance.",
        "scope": {
            "scope_type": "broad",
            "compatible_functional_groups": ["sp2_halide_present"],
            "incompatible_functional_groups": [],
        },
        "applies_if": {
            "all": ["sp2_halide_present", "sp2_boron_present"],
            "any": [],
        },
        "default_rule": {
            "id": "auto_default",
            "description": "Baseline conditions",
            "conditions": {
                "pd_source": "PdCl2(dtbpf)",
                "ligand": "dtbpf",
                "base": "K3PO4 (aq)",
                "solvent": "1,4-dioxane/H2O 4:1",
                "temperature_C": "90-100",
                "time_h": "1-2",
            },
        },
        "base_rules": [
            {
                "id": "auto_rule",
                "name": "Auto Rule",
                "description": "LLM-generated rule",
                "reactant_features": {"any": ["sp2_chloride_present"]},
                "conditions": {
                    "pd_source": "PdCl2(dtbpf)",
                    "ligand": "dtbpf",
                    "base": "K3PO4",
                    "solvent": "1,4-dioxane/H2O 4:1",
                    "temperature_C": "100",
                    "time_h": "2",
                },
                "priority": 0,
            }
        ],
        "modifiers": [
            {
                "id": "mod1",
                "when": ["symptom:hydrodehalogenation_observed"],
                "suggest": "Add TBAB 0.5-1.0 equiv",
                "rationale": "Suppress β-hydride pathways",
            }
        ],
    }

    monkeypatch.setattr(
        "chem_assistant.chemtools_wrapper._invoke_rule_builder_llm",
        lambda prompt: json.dumps(payload),
    )

    params = RuleBuilderAutoInput(
        family="Suzuki_Miyaura",
        metadata_id="auto_suzuki",
        metadata_name="Auto Suzuki Draft",
        metadata_version="v0",
        reference_reactions=[
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        ],
        protocol_text="HTE plates favored dtbpf/K3PO4 with aqueous co-solvent.",
    )

    result = run_rule_builder_autofill(params)
    assert result["success"] is True
    rule_db = result["rule_database"]
    assert rule_db["metadata"]["id"] == "auto_suzuki"
    assert rule_db["base_rules"], "Expected base rules to be populated"
    assert all("field" in issue for issue in result["issues"])


def test_rule_builder_autofill_bad_json(monkeypatch):
    monkeypatch.setattr(
        "chem_assistant.chemtools_wrapper._invoke_rule_builder_llm",
        lambda prompt: "not-json",
    )

    params = RuleBuilderAutoInput(
        family="Suzuki_Miyaura",
        metadata_id="auto_failure",
        metadata_name="Failure Draft",
        metadata_version="v0",
        reference_reactions=[
            "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        ],
        protocol_text="Placeholder text",
    )

    result = run_rule_builder_autofill(params)
    assert result["success"] is False
    assert "raw_response" in result
