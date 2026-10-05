"""Source mass/amount conflicts retain reports and assignment uncertainty."""

from dataclasses import asdict

import pytest

from condition_registry.quantity_audit import audit_source_quantities, quantity_consistency_rules
from condition_registry import condition_registry_definition_versions


def test_reported_methyl_iodide_discrepancy_is_not_silently_corrected() -> None:
    text = "Methyl iodide (162 g, 2.60 mol) was charged."
    result, = audit_source_quantities("CI", text, ("Methyl iodide",))
    assert result.status == "conflicting"
    assert result.reported_mass == 162
    assert result.reported_amount == 2.60
    assert result.expected_mass_g == pytest.approx(369.0414)
    assert text[result.source_start:result.source_end] == result.source_text
    assert "Conditional" in result.reason
    assert asdict(result) == asdict(audit_source_quantities("CI", text, ("Methyl iodide",))[0])


@pytest.mark.parametrize("text", ["Ethanol (46.1 g, 1 mol)", "Ethanol (1000 mmol, 46.1 g)",
                                      "Ethanol (46.1 mg, 1 mmol)", "Ethanol (46.1 µg, 1 μmol)"])
def test_rounded_consistent_quantities_and_si_units(text: str) -> None:
    result, = audit_source_quantities("CCO", text, ("Ethanol",))
    assert result.status == "consistent"


@pytest.mark.parametrize("text", ["Ethanol (1 mL, 1 mmol)", "Ethanol (46.1 g, 1 mol, 50%)",
                                      "Ethanol solution (46.1 g, 1 mol)",
                                      "Ethanol was charged. Water (46.1 g, 1 mol)"])
def test_unsupported_or_nonadjacent_reports_are_not_guessed(text: str) -> None:
    assert audit_source_quantities("CCO", text, ("Ethanol",)) == ()


@pytest.mark.parametrize("smiles", ["C1", "*", "[CH3]", "CCO.C"])
def test_invalid_or_ambiguous_material_retains_reports_without_accepting_conversion(smiles: str) -> None:
    result, = audit_source_quantities(smiles, "Material (46.1 g, 1 mol)", ("Material",))
    assert result.status == "not_assessed" and result.expected_mass_g is None


def test_full_salt_graph_is_not_stripped_and_aliases_do_not_duplicate_quantities() -> None:
    result, = audit_source_quantities("[NH4+].[Cl-]", "Ammonium chloride (53.49 g, 1 mol)",
                                     ("Ammonium chloride", "chloride"))
    assert result.status == "consistent"
    assert result.molecular_weight == pytest.approx(53.492)
    assert "quantity_consistency.v1.json" in condition_registry_definition_versions()
    assert quantity_consistency_rules()["relative_tolerance"] == 0.05


def test_conflicting_repeated_reports_are_both_retained() -> None:
    checks = audit_source_quantities("CCO", "Ethanol (46.1 g, 1 mol). Ethanol (10 g, 1 mol).", ("Ethanol",))
    assert [item.status for item in checks] == ["consistent", "conflicting"]


@pytest.mark.parametrize("text", ["Ethanol (1e309 g, 1 mol)", "Ethanol (1e308 g, 1e-300 mol)"])
def test_nonfinite_reports_preserve_literal_evidence_without_nonfinite_json_values(text: str) -> None:
    result, = audit_source_quantities("CCO", text, ("Ethanol",))
    assert result.status == "not_assessed" and result.relative_error is None
    assert result.source_text == text
    import json
    json.dumps(asdict(result), allow_nan=False)


@pytest.mark.parametrize("change", [{"relative_tolerance": 0}, {"relative_tolerance": True},
                                     {"mass_units_g": {}}, {"amount_units_mol": {"mol": -1}},
                                     {"schema_version": "unknown"}])
def test_quantity_definition_loader_rejects_invalid_policies(tmp_path, monkeypatch, change) -> None:
    import json
    import condition_registry.quantity_audit as module

    rules = {**quantity_consistency_rules(), **change}
    directory = tmp_path / "definitions"
    directory.mkdir()
    (directory / "quantity_consistency.v1.json").write_text(json.dumps(rules), encoding="utf-8")
    monkeypatch.setattr(module, "__file__", str(tmp_path / "quantity_audit.py"))
    module.quantity_consistency_rules.cache_clear()
    try:
        with pytest.raises(ValueError):
            module.quantity_consistency_rules()
    finally:
        module.quantity_consistency_rules.cache_clear()
