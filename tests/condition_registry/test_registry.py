from condition_registry import (
    get_registry,
    load_condition_vocabulary,
    resolve_condition_components,
    resolve_substance,
    resolve_substance_id,
    validate_registry,
)
from condition_registry.normalization import normalize_cas


def test_exact_cas_resolution_preserves_multiple_roles() -> None:
    result = resolve_substance(cas="14221-01-3")
    assert result.status == "resolved"
    assert result.match_kind == "exact_cas"
    assert result.substance.canonical_name == "Tetrakis(triphenylphosphine)palladium"
    assert "metal_catalyst" in {role.role_id for role in result.substance.roles}


def test_invalid_cas_is_not_silently_resolved() -> None:
    assert normalize_cas("7732-18-4") is None
    result = resolve_substance(cas="7732-18-4")
    assert result.status == "invalid_identifier"


def test_exact_substance_id_resolution_is_available_to_definitions() -> None:
    result = resolve_substance_id("cas:1536473-72-9")

    assert result.status == "resolved"
    assert result.match_kind == "exact_substance_id"
    assert result.substance is not None
    assert result.substance.canonical_name == "t-BuBrettPhos Palladacycle Gen. 3"


def test_condition_component_resolution_keeps_source_fields() -> None:
    resolved = resolve_condition_components(
        catalyst_cas=("14221-01-3",),
        reagent_cas=("584-08-7",),
        solvent_cas=("7732-18-5",),
    )
    assert set(resolved) == {"catalyst_cas", "reagent_cas", "solvent_cas"}
    assert all(items[0].status == "resolved" for items in resolved.values())


def test_registry_audit_reconciles_all_rows() -> None:
    report = validate_registry()
    assert report["total_rows"] == get_registry().size
    assert report["accepted_rows"] + report["issue_rows"] == report["total_rows"]
    assert report["issue_rows"] > 0


def test_condition_vocabulary_is_immutable_and_versioned() -> None:
    vocabulary = load_condition_vocabulary()

    assert "metal_catalyst" in vocabulary.role_ids
    assert "pd_zero_sources" in vocabulary.family_ids
    assert vocabulary.schema_version == "roles_families.v1"
