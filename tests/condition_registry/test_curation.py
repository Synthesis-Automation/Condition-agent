import json
from pathlib import Path

import pytest

from condition_registry import (
    CompoundAdditionError,
    CompoundAdditionRequest,
    CompoundAliasInput,
    ConditionRegistry,
    RoleCapability,
    SubstanceAliasAdditionRequest,
    add_compound,
    add_substance_aliases,
    update_compound,
    validate_registry,
)
from condition_registry import curation


def _water_record() -> dict:
    return {
        "id": "cas:7732-18-5",
        "name": "Water",
        "cas": "7732-18-5",
        "smiles": "O",
        "roles": ["solvent"],
    }


def _definition_path(tmp_path: Path) -> Path:
    path = tmp_path / "substances.v2.jsonl"
    path.write_text(json.dumps(_water_record()) + "\n", encoding="utf-8")
    return path


def _request(**changes) -> CompoundAdditionRequest:
    values = {
        "canonical_name": "Ethanol",
        "cas": "64-17-5",
        "smiles": "CCO",
        "abbreviation": "EtOH",
        "substance_kind": "liquid",
        "density": 0.789,
        "roles": (
            RoleCapability("solvent", tag="protic solvent", evidence="curator"),
            RoleCapability("additive", evidence="curator"),
            RoleCapability("reductant", evidence="curator"),
        ),
        "aliases": (
            CompoundAliasInput("systematic_name", "Ethyl alcohol", "en"),
            CompoundAliasInput("common_name", "Spirit of wine", "en"),
        ),
    }
    values.update(changes)
    return CompoundAdditionRequest(**values)


def test_add_compound_writes_one_record_with_nested_aliases_and_roles(tmp_path) -> None:
    path = _definition_path(tmp_path)

    result = add_compound(_request(), substances_path=path)

    assert result.substance.substance_id == "cas:64-17-5"
    assert result.formula == "C2H6O"
    assert result.molecular_weight == pytest.approx(46.069, abs=0.001)
    assert {role.role_id for role in result.substance.roles} == {
        "solvent",
        "additive",
        "reductant",
    }
    records = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    assert len(records) == 2
    assert len(records[1]["aliases"]) == 3
    assert records[1]["roles"] == ["solvent", "additive", "reductant"]
    assert set(records[1]) == {"id", "name", "cas", "smiles", "aliases", "roles"}
    registry = ConditionRegistry(substances_path=path)
    assert registry.resolve(name="Ethyl alcohol").substance == result.substance


def test_add_compound_rejects_duplicate_cas_without_writing(tmp_path) -> None:
    path = _definition_path(tmp_path)
    before = path.read_bytes()

    with pytest.raises(CompoundAdditionError, match="CAS_ALREADY_REGISTERED"):
        add_compound(
            _request(canonical_name="Duplicate water", cas="7732-18-5"),
            substances_path=path,
        )

    assert path.read_bytes() == before


def test_add_compound_rejects_structure_and_unknown_role(tmp_path) -> None:
    path = _definition_path(tmp_path)

    with pytest.raises(CompoundAdditionError) as caught:
        add_compound(
            _request(
                formula="C3H8O",
                roles=(RoleCapability("invented_role"),),
            ),
            substances_path=path,
        )

    assert "FORMULA_SMILES_MISMATCH:C3H8O:C2H6O" in caught.value.errors
    assert "UNKNOWN_ROLE:invented_role" in caught.value.errors


def test_add_compound_rolls_back_unified_definition_on_failure(tmp_path, monkeypatch) -> None:
    path = _definition_path(tmp_path)
    before = path.read_bytes()

    def fail_write(path, records):
        path.write_text("partial", encoding="utf-8")
        raise OSError("simulated unified write failure")

    monkeypatch.setattr(curation, "_write_jsonl_atomic", fail_write)
    with pytest.raises(OSError, match="simulated unified write failure"):
        add_compound(_request(), substances_path=path)

    assert path.read_bytes() == before


def test_shared_alias_is_explicitly_ambiguous_in_both_records(tmp_path) -> None:
    path = _definition_path(tmp_path)
    add_compound(
        _request(aliases=(CompoundAliasInput("common_name", "Shared spirit"),)),
        substances_path=path,
    )
    add_compound(
        _request(
            canonical_name="Methanol",
            cas="67-56-1",
            smiles="CO",
            abbreviation="MeOH",
            aliases=(
                CompoundAliasInput(
                    "common_name", "Shared spirit", allow_ambiguous=True
                ),
            ),
        ),
        substances_path=path,
    )

    result = ConditionRegistry(substances_path=path).resolve(name="Shared spirit")
    assert result.status == "ambiguous"
    assert result.candidates == ("cas:64-17-5", "cas:67-56-1")
    assert not validate_registry(substances_path=path)["has_errors"]


def test_update_compound_replaces_nested_aliases_and_roles(tmp_path) -> None:
    path = _definition_path(tmp_path)
    created = add_compound(_request(), substances_path=path)

    updated = update_compound(
        created.substance.substance_id,
        _request(
            density=0.79,
            roles=(RoleCapability("solvent"),),
            aliases=(CompoundAliasInput("common_name", "Alcohol"),),
        ),
        substances_path=path,
    )

    assert updated.substance.properties == {}
    assert updated.substance.provenance == {}
    assert updated.substance.aliases == ("EtOH", "Alcohol")
    registry = ConditionRegistry(substances_path=path)
    assert registry.resolve(name="Spirit of wine").status == "unresolved"


def test_update_rejects_cas_change_without_writing(tmp_path) -> None:
    path = _definition_path(tmp_path)
    before = path.read_bytes()

    with pytest.raises(CompoundAdditionError, match="CAS_CHANGE_NOT_ALLOWED"):
        update_compound(
            "cas:7732-18-5",
            CompoundAdditionRequest(
                canonical_name="Water", cas="64-17-5", smiles="O"
            ),
            substances_path=path,
        )

    assert path.read_bytes() == before


def test_add_substance_aliases_updates_target_record_atomically(tmp_path) -> None:
    path = _definition_path(tmp_path)

    result = add_substance_aliases(
        (
            SubstanceAliasAdditionRequest(
                substance_id="cas:7732-18-5",
                identifier_type="legacy_name",
                value="Aqua",
                language="en",
            ),
        ),
        substances_path=path,
    )

    assert tuple(item.value for item in result.added) == ("Aqua",)
    resolved = ConditionRegistry(substances_path=path).resolve(name="Aqua")
    assert resolved.substance is not None
    assert resolved.substance.substance_id == "cas:7732-18-5"
