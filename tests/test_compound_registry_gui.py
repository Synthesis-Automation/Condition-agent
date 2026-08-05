from app import compound_registry_gui as gui
from condition_registry import (
    CompoundAdditionError,
    CompoundAdditionResult,
    ResolutionResult,
    Substance,
)


def _set_alias_value(
    window: gui.CompoundRegistryWindow,
    row: int,
    value: str,
) -> None:
    value_edit = window.alias_table.cellWidget(row, 1)
    value_edit.setText(value)


def test_compound_registry_window_exposes_curated_fields(qtbot) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)

    assert window.windowTitle() == "Condition Registry — Add Compound"
    assert window.cas_edit.objectName() == "casNumber"
    assert window.name_edit.objectName() == "canonicalName"
    assert window.source_edit.objectName() == "provenanceSource"
    assert window.alias_table.rowCount() == 1
    assert window.primary_role_combo.currentData() == "other_reagent"
    assert window.primary_family_combo.findData("misc_general") >= 0
    assert window.status_box.isReadOnly()


def test_compound_form_builds_typed_request_and_filters_families(qtbot) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    window.cas_edit.setText("64-17-5")
    window.name_edit.setText("Ethanol")
    window.abbreviation_edit.setText("EtOH")
    window.smiles_edit.setText("CCO")
    window.source_edit.setText("chemist:test")
    window.density_edit.setText("0.789")
    role_index = window.primary_role_combo.findData("solvent")
    window.primary_role_combo.setCurrentIndex(role_index)
    family_index = window.primary_family_combo.findData("alcohols_primary")
    window.primary_family_combo.setCurrentIndex(family_index)
    _set_alias_value(window, 0, "Ethyl alcohol")

    request = window.compound_request()

    assert request.cas == "64-17-5"
    assert request.canonical_name == "Ethanol"
    assert request.density == 0.789
    assert request.roles[0].role_id == "solvent"
    assert request.roles[0].family_id == "alcohols_primary"
    assert request.aliases[0].identifier_type == "common_name"
    assert request.aliases[0].value == "Ethyl alcohol"


def test_save_compound_calls_registry_service_and_reports_result(
    qtbot,
    monkeypatch,
) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    window.cas_edit.setText("64-17-5")
    window.name_edit.setText("Ethanol")
    window.smiles_edit.setText("CCO")
    window.source_edit.setText("chemist:test")
    _set_alias_value(window, 0, "Ethyl alcohol")
    captured = []
    substance = Substance(
        substance_id="cas:64-17-5",
        canonical_name="Ethanol",
        cas="64-17-5",
        smiles="CCO",
    )

    def fake_add(request):
        captured.append(request)
        return CompoundAdditionResult(
            substance=substance,
            canonical_smiles="CCO",
            formula="C2H6O",
            molecular_weight=46.069,
            alias_count=1,
            additions_path="additions.csv",
            identifiers_path="identifiers.csv",
        )

    monkeypatch.setattr(gui, "add_compound", fake_add)

    window.save_compound()

    assert len(captured) == 1
    assert captured[0].source == "chemist:test"
    assert "Added Ethanol as cas:64-17-5" in window.status_box.toPlainText()
    assert window.cas_status.text() == "Added"
    assert window.save_button.isEnabled()


def test_save_compound_displays_all_validation_errors(qtbot, monkeypatch) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)

    def reject(_request):
        raise CompoundAdditionError(("INVALID_CAS", "MISSING_SOURCE"))

    monkeypatch.setattr(gui, "add_compound", reject)

    window.save_compound()

    status = window.status_box.toPlainText()
    assert "INVALID_CAS" in status
    assert "MISSING_SOURCE" in status
    assert window.cas_status.text() == "Validation failed"


def test_check_cas_reports_existing_identity(qtbot, monkeypatch) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    window.cas_edit.setText("7732-18-5")
    water = Substance(
        substance_id="cas:7732-18-5",
        canonical_name="Water",
        cas="7732-18-5",
        smiles="O",
    )
    monkeypatch.setattr(
        gui,
        "resolve_substance",
        lambda **_query: ResolutionResult(
            query="7732-18-5",
            status="resolved",
            substance=water,
        ),
    )

    window.check_cas()

    assert window.cas_status.text() == "Exists: Water"
