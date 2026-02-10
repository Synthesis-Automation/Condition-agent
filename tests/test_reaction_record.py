from chemtools.smiles import normalize_reaction as normalize_reaction_legacy
from chemtools.featurizers.analysis.reaction_record import ReactionRecord
from chemtools.featurizers.analysis.smiles import normalize_reaction


def test_reaction_record_round_trip_payload() -> None:
    payload = normalize_reaction("CCO.CBr>>CCOC")
    record = ReactionRecord.from_payload(payload)

    assert record.reactant_smiles
    assert record.product_smiles
    assert record.to_payload() == payload


def test_normalize_reaction_handles_nonstandard_triple_arrow_format() -> None:
    payload = normalize_reaction("CCO>>O>>CCOC")
    record = ReactionRecord.from_payload(payload)

    # Non-standard reactants>>reagent>>products is normalized into reactants + reagent.
    assert len(record.reactants) == 2
    assert len(record.agents) == 0
    assert len(record.products) == 1
    assert payload["normalized"] == record.normalized


def test_normalize_reaction_contract_keys_and_compatibility() -> None:
    reaction = "CCO.CBr>>CCOC"
    payload = normalize_reaction(reaction)
    legacy_payload = normalize_reaction_legacy(reaction)

    assert payload == legacy_payload
    assert set(payload.keys()) == {
        "input",
        "reactants",
        "agents",
        "products",
        "normalized",
        "errors",
    }
    assert isinstance(payload["reactants"], list)
    assert isinstance(payload["agents"], list)
    assert isinstance(payload["products"], list)
    assert isinstance(payload["errors"], list)
    assert isinstance(payload["normalized"], str)
    assert payload["normalized"] == ReactionRecord.from_payload(payload).normalized


def test_normalize_reaction_contract_reports_invalid_component_errors() -> None:
    payload = normalize_reaction("CCO.bad>>CCO")

    assert "bad" in payload["errors"]
    bad_entries = [
        component
        for component in payload["reactants"]
        if component.get("input") == "bad"
    ]
    assert bad_entries
    assert bad_entries[0].get("error") == "INVALID_SMILES"
