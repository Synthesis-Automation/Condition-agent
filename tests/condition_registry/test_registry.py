from condition_registry import (
    ConditionRegistry,
    Substance,
    SubstanceIdentifier,
    get_registry,
    load_condition_vocabulary,
    resolve_condition_components,
    resolve_identifier,
    resolve_substance,
    resolve_substance_id,
    validate_registry,
)
from condition_registry.normalization import normalize_cas


def _identifier(
    identifier_id: str,
    substance_id: str,
    identifier_type: str,
    value: str,
    *,
    status: str = "active",
) -> SubstanceIdentifier:
    return SubstanceIdentifier(
        identifier_id=identifier_id,
        substance_id=substance_id,
        identifier_type=identifier_type,
        value=value,
        source="test_curator",
        status=status,
    )


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
    assert report["identifier_total_rows"] >= 2
    assert report["identifier_issue_rows"] == 0


def test_condition_vocabulary_is_immutable_and_versioned() -> None:
    vocabulary = load_condition_vocabulary()

    assert "metal_catalyst" in vocabulary.role_ids
    assert "pd_zero_sources" in vocabulary.family_ids
    assert vocabulary.schema_version == "roles_families.v1"


def test_substance_supports_unbounded_typed_aliases_with_provenance() -> None:
    result = resolve_identifier(
        "Dipotassium carbonate",
        identifier_type="systematic_name",
    )

    assert result.status == "resolved"
    assert result.substance is not None
    assert result.substance.substance_id == "cas:584-08-7"
    assert {"K2CO3", "Dipotassium carbonate"} <= set(
        result.substance.aliases
    )
    assert result.matched_identifier is not None
    assert result.matched_identifier.source == "condition_registry_curated"
    assert result.matched_identifier.identifier_type == "systematic_name"


def test_shared_alias_returns_all_candidates_instead_of_guessing() -> None:
    first_id = "sub:test:first"
    second_id = "sub:test:second"
    first = Substance(
        substance_id=first_id,
        canonical_name="First compound",
        cas=None,
        smiles="C",
        identifiers=(
            _identifier(
                "id:first:name", first_id, "canonical_name", "First compound"
            ),
            _identifier("id:first:shared", first_id, "common_name", "Shared alias"),
        ),
    )
    second = Substance(
        substance_id=second_id,
        canonical_name="Second compound",
        cas=None,
        smiles="CC",
        identifiers=(
            _identifier(
                "id:second:name", second_id, "canonical_name", "Second compound"
            ),
            _identifier("id:second:shared", second_id, "common_name", "Shared alias"),
        ),
    )

    result = ConditionRegistry(substances=(first, second)).resolve_identifier(
        "shared alias",
        identifier_type="common_name",
    )

    assert result.status == "ambiguous"
    assert result.candidates == (first_id, second_id)
    assert result.substance is None
    assert result.warnings == ("AMBIGUOUS_IDENTIFIER",)


def test_deprecated_alias_is_retained_but_not_used_for_resolution() -> None:
    substance_id = "sub:test:deprecated"
    substance = Substance(
        substance_id=substance_id,
        canonical_name="Current name",
        cas=None,
        smiles="N",
        identifiers=(
            _identifier("id:current", substance_id, "canonical_name", "Current name"),
            _identifier(
                "id:deprecated",
                substance_id,
                "legacy_name",
                "Retired name",
                status="deprecated",
            ),
        ),
    )
    registry = ConditionRegistry(substances=(substance,))

    assert registry.resolve(name="Current name").status == "resolved"
    assert registry.resolve(name="Retired name").status == "unresolved"
    assert substance.identifiers[1].status == "deprecated"
    assert "Retired name" not in substance.aliases


def test_formulation_words_are_not_erased_by_typed_name_resolution() -> None:
    solution_id = "sub:test:solution"
    dry_id = "sub:test:dry"
    solution = Substance(
        substance_id=solution_id,
        canonical_name="Example solution",
        cas=None,
        smiles=None,
        identifiers=(
            _identifier(
                "id:solution",
                solution_id,
                "canonical_name",
                "Example reagent aqueous solution",
            ),
        ),
    )
    dry = Substance(
        substance_id=dry_id,
        canonical_name="Dry example",
        cas=None,
        smiles=None,
        identifiers=(
            _identifier(
                "id:dry",
                dry_id,
                "canonical_name",
                "Example reagent anhydrous",
            ),
        ),
    )
    registry = ConditionRegistry(substances=(solution, dry))

    assert registry.resolve(name="Example reagent aqueous solution").substance == solution
    assert registry.resolve(name="Example reagent anhydrous").substance == dry


def test_identifier_validation_rejects_unknown_substance_and_missing_source(
    tmp_path,
) -> None:
    identifiers_path = tmp_path / "identifiers.csv"
    identifiers_path.write_text(
        "identifier_id,substance_id,identifier_type,value,language,"
        "is_preferred,source,confidence,status,normalization_profile,"
        "allow_ambiguous\n"
        "sid:unknown,sub:missing,common_name,Unknown alias,en,false,,1.0,"
        "active,chemical_name_v1,false\n",
        encoding="utf-8",
    )

    report = validate_registry(identifiers_path=identifiers_path)

    assert report["identifier_issue_rows"] == 1
    assert {
        "UNKNOWN_SUBSTANCE_ID:sub:missing",
        "MISSING_IDENTIFIER_SOURCE",
    } <= set(report["identifier_issues"][0]["issues"])


def test_identifier_validation_requires_declared_shared_aliases(tmp_path) -> None:
    identifiers_path = tmp_path / "identifiers.csv"
    header = (
        "identifier_id,substance_id,identifier_type,value,language,"
        "is_preferred,source,confidence,status,normalization_profile,"
        "allow_ambiguous\n"
    )
    identifiers_path.write_text(
        header
        + "sid:shared:one,cas:584-08-7,common_name,Shared curated alias,en,"
        "false,test_curator,1.0,active,chemical_name_v1,false\n"
        + "sid:shared:two,cas:7732-18-5,common_name,Shared curated alias,en,"
        "false,test_curator,1.0,active,chemical_name_v1,false\n",
        encoding="utf-8",
    )

    report = validate_registry(identifiers_path=identifiers_path)

    assert report["identifier_issue_rows"] == 2
    assert all(
        "UNDECLARED_AMBIGUOUS_IDENTIFIER" in item["issues"]
        for item in report["identifier_issues"]
    )
