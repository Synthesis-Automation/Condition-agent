import csv
from pathlib import Path

import pytest

from condition_registry import (
    CompoundAdditionError,
    CompoundAdditionRequest,
    CompoundAliasInput,
    ConditionRegistry,
    RoleAssignment,
    add_compound,
    update_compound,
    validate_registry,
)
from condition_registry import curation


LEGACY_FIELDNAMES = (
    "name",
    "abbreviation",
    "cas",
    "smiles",
    "formula",
    "type",
    "density",
    "mw",
    "bp",
    "mp",
    "volatile",
    "viscose",
    "role_1",
    "family_1",
    "tag_1",
    "role_2",
    "family_2",
    "tag_2",
)


def _write_csv(
    path: Path,
    fieldnames: tuple[str, ...],
    rows: tuple[dict[str, str], ...] = (),
) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _definition_paths(tmp_path: Path) -> tuple[Path, Path, Path]:
    substances = tmp_path / "substances.v1.csv"
    additions = tmp_path / "substance_additions.v1.csv"
    identifiers = tmp_path / "substance_identifiers.v1.csv"
    _write_csv(
        substances,
        LEGACY_FIELDNAMES,
        (
            {
                "name": "Water",
                "abbreviation": "H2O",
                "cas": "7732-18-5",
                "smiles": "O",
                "formula": "H2O",
                "type": "liquid",
                "role_1": "solvent",
                "family_1": "water",
            },
        ),
    )
    _write_csv(additions, curation.ADDITION_FIELDNAMES)
    _write_csv(identifiers, curation.IDENTIFIER_FIELDNAMES)
    return substances, additions, identifiers


def _request(**changes) -> CompoundAdditionRequest:
    values = {
        "canonical_name": "Ethanol",
        "cas": "64-17-5",
        "source": "chemist:test",
        "smiles": "CCO",
        "abbreviation": "EtOH",
        "substance_kind": "liquid",
        "density": 0.789,
        "boiling_point_c": 78.37,
        "melting_point_c": -114.1,
        "roles": (
            RoleAssignment("solvent", "alcohols_primary", "protic solvent"),
        ),
        "aliases": (
            CompoundAliasInput("systematic_name", "Ethyl alcohol", "en"),
            CompoundAliasInput("common_name", "Spirit of wine", "en"),
        ),
        "curator_notes": "Added through the curation workflow.",
    }
    values.update(changes)
    return CompoundAdditionRequest(**values)


def test_add_compound_writes_identity_and_unbounded_alias_rows(tmp_path) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)

    result = add_compound(
        _request(),
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )

    assert result.substance.substance_id == "cas:64-17-5"
    assert result.canonical_smiles == "CCO"
    assert result.formula == "C2H6O"
    assert result.molecular_weight == pytest.approx(46.069, abs=0.001)
    assert result.alias_count == 2
    registry = ConditionRegistry(
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    assert registry.resolve(cas="64-17-5").substance == result.substance
    alias = registry.resolve_identifier(
        "Ethyl alcohol", identifier_type="systematic_name"
    )
    assert alias.status == "resolved"
    assert alias.matched_identifier is not None
    assert alias.matched_identifier.source == "chemist:test"
    with additions.open(encoding="utf-8-sig", newline="") as handle:
        added = list(csv.DictReader(handle))
    assert added[0]["source"] == "chemist:test"
    assert added[0]["role_1"] == "solvent"
    assert added[0]["family_1"] == "alcohols_primary"


def test_add_compound_rejects_duplicate_cas_without_writing(tmp_path) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)
    before = additions.read_bytes()

    with pytest.raises(CompoundAdditionError, match="CAS_ALREADY_REGISTERED"):
        add_compound(
            _request(canonical_name="Duplicate water", cas="7732-18-5"),
            substances_path=substances,
            additions_path=additions,
            identifiers_path=identifiers,
        )

    assert additions.read_bytes() == before


def test_add_compound_rejects_structure_and_taxonomy_conflicts(tmp_path) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)

    with pytest.raises(CompoundAdditionError) as caught:
        add_compound(
            _request(
                formula="C3H8O",
                roles=(RoleAssignment("base", "water"),),
            ),
            substances_path=substances,
            additions_path=additions,
            identifiers_path=identifiers,
        )

    assert "FORMULA_SMILES_MISMATCH:C3H8O:C2H6O" in caught.value.errors
    assert "ROLE_FAMILY_MISMATCH:base:water" in caught.value.errors
    assert additions.read_text(encoding="utf-8-sig").count("\n") == 1


def test_add_compound_rolls_back_both_definitions_on_write_failure(
    tmp_path,
    monkeypatch,
) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)
    additions_before = additions.read_bytes()
    identifiers_before = identifiers.read_bytes()
    real_write = curation._write_csv_atomic
    calls = 0

    def fail_second_write(path, fieldnames, rows):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated identifier write failure")
        real_write(path, fieldnames, rows)

    monkeypatch.setattr(curation, "_write_csv_atomic", fail_second_write)

    with pytest.raises(OSError, match="simulated identifier write failure"):
        add_compound(
            _request(),
            substances_path=substances,
            additions_path=additions,
            identifiers_path=identifiers,
        )

    assert additions.read_bytes() == additions_before
    assert identifiers.read_bytes() == identifiers_before


def test_add_compound_can_explicitly_share_an_existing_curated_alias(
    tmp_path,
) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)
    first = _request(
        aliases=(CompoundAliasInput("common_name", "Shared spirit", "en"),)
    )
    add_compound(
        first,
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    second = _request(
        canonical_name="Methanol",
        cas="67-56-1",
        smiles="CO",
        abbreviation="MeOH",
        aliases=(
            CompoundAliasInput(
                "common_name",
                "Shared spirit",
                "en",
                allow_ambiguous=True,
            ),
        ),
    )

    add_compound(
        second,
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )

    registry = ConditionRegistry(
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    result = registry.resolve(name="Shared spirit")
    assert result.status == "ambiguous"
    assert result.candidates == ("cas:64-17-5", "cas:67-56-1")
    with identifiers.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["allow_ambiguous"] for row in rows} == {"true"}
    report = validate_registry(
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    assert report["issue_rows"] == 0
    assert report["identifier_issue_rows"] == 0


def test_update_compound_replaces_addition_and_its_aliases(tmp_path) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)
    created = add_compound(
        _request(),
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    updated_request = _request(
        density=0.79,
        aliases=(
            CompoundAliasInput("common_name", "Alcohol", "en"),
        ),
        curator_notes="Reviewed and updated.",
        source="chemist:update",
    )

    updated = update_compound(
        created.substance.substance_id,
        updated_request,
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )

    assert updated.substance.substance_id == "cas:64-17-5"
    assert updated.substance.properties["density"] == "0.79"
    assert updated.substance.properties["source"] == "chemist:update"
    assert updated.substance.aliases == ("EtOH", "Alcohol")
    registry = ConditionRegistry(
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    assert registry.resolve(name="Spirit of wine").status == "unresolved"
    assert registry.resolve(name="Alcohol").substance == updated.substance


def test_update_compound_migrates_legacy_row_to_curated_additions(
    tmp_path,
) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)
    request = CompoundAdditionRequest(
        canonical_name="Water",
        cas="7732-18-5",
        source="chemist:legacy-review",
        smiles="O",
        abbreviation="H2O",
        substance_kind="liquid",
        density=0.997,
        roles=(RoleAssignment("solvent", "water"),),
        aliases=(CompoundAliasInput("systematic_name", "Oxidane", "en"),),
        curator_notes="Migrated through the edit workflow.",
    )

    result = update_compound(
        "cas:7732-18-5",
        request,
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )

    with substances.open(encoding="utf-8-sig", newline="") as handle:
        assert list(csv.DictReader(handle)) == []
    with additions.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["substance_id"] == "cas:7732-18-5"
    assert rows[0]["source"] == "chemist:legacy-review"
    assert result.substance.canonical_name == "Water"
    assert result.substance.aliases == ("H2O", "Oxidane")
    report = validate_registry(
        substances_path=substances,
        additions_path=additions,
        identifiers_path=identifiers,
    )
    assert report["issue_rows"] == 0
    assert report["identifier_issue_rows"] == 0


def test_update_compound_rejects_cas_change_without_writing(tmp_path) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)
    before = (substances.read_bytes(), additions.read_bytes(), identifiers.read_bytes())

    with pytest.raises(CompoundAdditionError, match="CAS_CHANGE_NOT_ALLOWED"):
        update_compound(
            "cas:7732-18-5",
            CompoundAdditionRequest(
                canonical_name="Water",
                cas="64-17-5",
                source="chemist:test",
                smiles="O",
            ),
            substances_path=substances,
            additions_path=additions,
            identifiers_path=identifiers,
        )

    assert (
        substances.read_bytes(),
        additions.read_bytes(),
        identifiers.read_bytes(),
    ) == before


def test_update_compound_rolls_back_all_definitions_on_failure(
    tmp_path,
    monkeypatch,
) -> None:
    substances, additions, identifiers = _definition_paths(tmp_path)
    before = (substances.read_bytes(), additions.read_bytes(), identifiers.read_bytes())
    real_write = curation._write_csv_atomic
    calls = 0

    def fail_identifier_write(path, fieldnames, rows):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise OSError("simulated update failure")
        real_write(path, fieldnames, rows)

    monkeypatch.setattr(curation, "_write_csv_atomic", fail_identifier_write)

    with pytest.raises(OSError, match="simulated update failure"):
        update_compound(
            "cas:7732-18-5",
            CompoundAdditionRequest(
                canonical_name="Water",
                cas="7732-18-5",
                source="chemist:test",
                smiles="O",
            ),
            substances_path=substances,
            additions_path=additions,
            identifiers_path=identifiers,
        )

    assert (
        substances.read_bytes(),
        additions.read_bytes(),
        identifiers.read_bytes(),
    ) == before
