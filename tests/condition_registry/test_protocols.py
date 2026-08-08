from condition_registry import (
    ConditionComponentInput,
    ProtocolReactionMaterialInput,
    build_resolved_recipe_from_inputs,
    build_synthesis_protocol_draft,
)


def test_protocol_serializes_condition_cas_and_observed_setpoints() -> None:
    recipe = build_resolved_recipe_from_inputs(
        (
            ConditionComponentInput(
                "K2CO3",
                source_field="base",
                identifier_type="name",
                source_role_hint="base",
                amount=2.0,
                amount_unit="equiv",
            ),
        ),
        temperature_c=80.0,
        time_h=12.0,
        atmosphere="nitrogen",
    )
    protocol = build_synthesis_protocol_draft(
        recipe,
        reaction_smiles="CCBr.N>>CCN",
        reaction_materials=(
            ProtocolReactionMaterialInput("reactant", 0, "CCBr"),
            ProtocolReactionMaterialInput("reactant", 1, "N"),
            ProtocolReactionMaterialInput("product", 0, "CCN"),
        ),
    )
    payload = protocol.to_dict()

    condition = next(
        item for item in payload["materials"] if item["category"] == "condition"
    )
    assert condition["substance_id"] == "cas:584-08-7"
    assert condition["cas"] == "584-08-7"
    assert condition["role"] == "base"
    assert payload["operating_conditions"] == {
        "temperature_c": 80.0,
        "time_h": 12.0,
        "concentration_m": None,
        "atmosphere": "nitrogen",
    }
    assert payload["operations"][0]["operation_type"] == "maintain_conditions"
    assert payload["execution_readiness"] == "review_required"
    assert "ordered_operations" in payload["missing_required_fields"]


def test_protocol_id_is_deterministic_and_missing_data_is_explicit() -> None:
    recipe = build_resolved_recipe_from_inputs(
        (
            ConditionComponentInput(
                "water",
                source_field="solvent",
                identifier_type="name",
                source_role_hint="solvent",
            ),
        )
    )

    first = build_synthesis_protocol_draft(recipe)
    second = build_synthesis_protocol_draft(recipe)

    assert first.protocol_id == second.protocol_id
    assert first.protocol_id.startswith("CP1:")
    assert "reaction_inputs" in first.missing_required_fields
    assert "materials.condition_001.amount" in first.missing_required_fields
    assert "operating_conditions.temperature_c" in first.missing_required_fields
    assert "operating_conditions.time_h" in first.missing_required_fields


def test_protocol_enriches_cas_for_compatible_index_recipe() -> None:
    protocol = build_synthesis_protocol_draft(
        {
            "recipe_id": "RCR2:legacy-compatible",
            "recipe_core_id": "RCORE2:legacy-compatible",
            "bases": (
                {
                    "identity_status": "resolved",
                    "substance_id": "cas:584-08-7",
                    "canonical_name": "Potassium carbonate",
                    "primary_role": "base",
                    "role_status": "assigned",
                    "source_field": "reagent_cas",
                },
            ),
        }
    )

    assert protocol.materials[0].cas == "584-08-7"
