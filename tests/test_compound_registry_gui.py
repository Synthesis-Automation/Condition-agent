from app import compound_registry_gui as gui
from cas_tools import CompoundLookupResult
from condition_registry import (
    CompoundAdditionError,
    CompoundAdditionResult,
    ResolutionResult,
    RoleAssignment,
    Substance,
    SubstanceIdentifier,
)


def _set_alias_value(
    window: gui.CompoundRegistryWindow,
    row: int,
    value: str,
) -> None:
    value_edit = window.alias_table.cellWidget(row, 1)
    value_edit.setText(value)


def _lookup_result() -> CompoundLookupResult:
    return CompoundLookupResult(
        cas="64-17-5",
        status="resolved",
        canonical_name="Ethanol",
        iupac_name="ethanol",
        abbreviation="EtOH",
        smiles="CCO",
        formula="C2H6O",
        molecular_weight=46.07,
        inchi_key="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        pubchem_cid=702,
        density=0.7893,
        boiling_point_c=78.2,
        melting_point_c=-114.1,
        substance_kind="liquid",
        synonyms=("ethyl alcohol", "grain alcohol"),
        source_ids=("pubchem:cid:702",),
        source_urls=("https://pubchem.ncbi.nlm.nih.gov/compound/702",),
    )


def test_compound_registry_window_exposes_curated_fields(qtbot) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)

    assert window.windowTitle() == "Condition Registry — Add Compound"
    assert window.isMaximized()
    assert window.cas_edit.objectName() == "casNumber"
    assert window.name_edit.objectName() == "canonicalName"
    assert window.source_edit.objectName() == "provenanceSource"
    assert window.lookup_button.objectName() == "autoFillCompoundButton"
    assert window.alias_table.rowCount() == 1
    assert window.primary_role_combo.currentData() == "other_reagent"
    assert window.primary_family_combo.findData("misc_general") >= 0
    assert window.status_box.isReadOnly()


def test_identity_and_role_capabilities_use_two_columns(qtbot) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)

    columns = window.findChild(
        gui.QtWidgets.QHBoxLayout,
        "identityAndRolesColumns",
    )

    assert columns is not None
    assert columns.count() == 2
    assert columns.itemAt(0).widget().objectName() == "identityGroup"
    assert (
        columns.itemAt(1).widget().objectName()
        == "conditionRoleCapabilitiesGroup"
    )


def test_aliases_and_provenance_use_two_columns(qtbot) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)

    columns = window.findChild(
        gui.QtWidgets.QHBoxLayout,
        "aliasesAndProvenanceColumns",
    )

    assert columns is not None
    assert columns.count() == 2
    assert columns.itemAt(0).widget().objectName() == "additionalAliasesGroup"
    assert columns.itemAt(1).widget().objectName() == "provenanceGroup"


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


def test_check_cas_loads_existing_identity_for_editing(qtbot, monkeypatch) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    window.cas_edit.setText("7732-18-5")
    water = Substance(
        substance_id="cas:7732-18-5",
        canonical_name="Water",
        cas="7732-18-5",
        smiles="O",
        roles=(RoleAssignment("solvent", "water", "aqueous medium"),),
        properties={
            "formula": "H2O",
            "type": "liquid",
            "density": "0.997",
            "mw": "18.015",
            "source": "chemist:existing",
            "curator_notes": "Previously reviewed.",
        },
        identifiers=(
            SubstanceIdentifier(
                identifier_id="legacy:water:abbreviation",
                substance_id="cas:7732-18-5",
                identifier_type="abbreviation",
                value="H2O",
            ),
            SubstanceIdentifier(
                identifier_id="sid:water:oxidane",
                substance_id="cas:7732-18-5",
                identifier_type="systematic_name",
                value="Oxidane",
                language="en",
                is_preferred=True,
                source="chemist:existing",
                normalization_profile="chemical_name_v1",
            ),
        ),
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

    assert window.cas_status.text() == "Loaded for editing"
    assert window.editing_substance_id == "cas:7732-18-5"
    assert window.cas_edit.isReadOnly()
    assert window.name_edit.text() == "Water"
    assert window.abbreviation_edit.text() == "H2O"
    assert window.formula_edit.text() == "H2O"
    assert window.density_edit.text() == "0.997"
    assert window.primary_role_combo.currentData() == "solvent"
    assert window.primary_family_combo.currentData() == "water"
    assert window.source_edit.text() == "chemist:existing"
    assert window.alias_inputs() == (
        gui.CompoundAliasInput(
            identifier_type="systematic_name",
            value="Oxidane",
            language="en",
            is_preferred=True,
        ),
    )
    assert window.save_button.text() == "Validate and Save Changes"


def test_save_existing_compound_calls_update_service(qtbot, monkeypatch) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    substance = Substance(
        substance_id="cas:7732-18-5",
        canonical_name="Water",
        cas="7732-18-5",
        smiles="O",
        properties={"source": "chemist:existing"},
    )
    window.load_existing_substance(substance)
    captured = []

    def fake_update(substance_id, request):
        captured.append((substance_id, request))
        return CompoundAdditionResult(
            substance=substance,
            canonical_smiles="O",
            formula="H2O",
            molecular_weight=18.015,
            alias_count=0,
            additions_path="additions.csv",
            identifiers_path="identifiers.csv",
        )

    monkeypatch.setattr(gui, "update_compound", fake_update)

    window.save_compound()

    assert captured[0][0] == "cas:7732-18-5"
    assert captured[0][1].canonical_name == "Water"
    assert window.cas_status.text() == "Updated"
    assert "Updated Water as cas:7732-18-5" in window.status_box.toPlainText()


def test_lookup_worker_forwards_typed_result(monkeypatch) -> None:
    expected = _lookup_result()
    monkeypatch.setattr(gui, "lookup_compound_by_cas", lambda cas: expected)
    worker = gui.CompoundLookupWorker("64-17-5")
    received = []
    worker.finished.connect(
        lambda result, error: received.append((result, error))
    )

    worker.run()

    assert received == [(expected, "")]


def test_apply_web_lookup_fills_supported_fields_and_preserves_role_review(
    qtbot,
) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    window.cas_edit.setText("64-17-5")

    window.apply_lookup_result(_lookup_result())

    assert window.name_edit.text() == "Ethanol"
    assert window.abbreviation_edit.text() == "EtOH"
    assert window.smiles_edit.text() == "CCO"
    assert window.formula_edit.text() == "C2H6O"
    assert window.molecular_weight_edit.text() == "46.07"
    assert window.density_edit.text() == "0.7893"
    assert window.boiling_point_edit.text() == "78.2"
    assert window.melting_point_edit.text() == "-114.1"
    assert window.kind_combo.currentText() == "liquid"
    assert window.source_edit.text() == "pubchem:cid:702"
    assert window.primary_role_combo.currentData() == "other_reagent"
    aliases = {(item.identifier_type, item.value) for item in window.alias_inputs()}
    assert aliases == {
        ("common_name", "ethyl alcohol"),
        ("common_name", "grain alcohol"),
        ("inchi_key", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"),
        ("database_id", "pubchem:cid:702"),
    }
    assert "Roles and families require curator review" in (
        window.status_box.toPlainText()
    )
    assert "https://pubchem.ncbi.nlm.nih.gov/compound/702" in (
        window.notes_edit.toPlainText()
    )


def test_invalid_cas_does_not_start_background_lookup(qtbot) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    window.cas_edit.setText("64-17-4")

    window.start_web_lookup()

    assert window.lookup_thread is None
    assert window.cas_status.text() == "Invalid CAS"
    assert "checksum is invalid" in window.status_box.toPlainText()


def test_auto_fill_button_runs_lookup_on_background_thread(
    qtbot,
    monkeypatch,
) -> None:
    window = gui.CompoundRegistryWindow()
    qtbot.addWidget(window)
    window.cas_edit.setText("64-17-5")
    monkeypatch.setattr(
        gui,
        "lookup_compound_by_cas",
        lambda _cas: _lookup_result(),
    )

    window.start_web_lookup()

    qtbot.waitUntil(lambda: window.last_lookup_result is not None, timeout=3000)
    qtbot.waitUntil(lambda: window.lookup_thread is None, timeout=3000)
    assert window.name_edit.text() == "Ethanol"
    assert window.lookup_button.isEnabled()
    assert window.save_button.isEnabled()
