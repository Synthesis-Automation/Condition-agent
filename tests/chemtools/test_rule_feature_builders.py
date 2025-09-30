from __future__ import annotations

from chemtools import recommend


def _mock_role_pack(acid: str, amine: str) -> dict:
    return {
        "assignments": {
            "electrophile": {"index": 0, "role": None, "matched_via_mask": False},
            "nucleophile": {"index": 1, "role": "amine", "matched_via_mask": True},
        },
        "reactants": [
            {"smiles": acid, "masks": {"amine": 0}},
            {"smiles": amine, "masks": {"amine": 1}},
        ],
    }


def test_amide_rule_feature_builder_enriches_payload() -> None:
    acid = "O=C(O)c1ccccc1"
    amine = "NCc1ccccc1"
    features = recommend.feat_molecular.featurize(acid, amine)
    role_pack = _mock_role_pack(acid, amine)

    payload = recommend._compose_rule_features(
        "Amide_Formation",
        features,
        role_pack,
        reactants=[acid, amine],
        detection={"rule_family": "Amide_Formation"},
    )

    assert payload["family"] == "Amide_Formation"
    assert payload["reaction_family"] == "amide_formation"
    assert payload["electrophile"]["class"] == "carboxylic acid"
    assert "carboxylic acid" in payload["substrate_class"].lower()
    assert payload.get("category")
    assert payload.get("water_management")
    assignments = payload.get("role_assignments") or {}
    assert assignments.get("electrophile", {}).get("smiles") == acid
    assert assignments.get("nucleophile", {}).get("smiles") == amine


def test_default_rule_feature_builder_remains_stable() -> None:
    electrophile = "Brc1ccccc1"
    nucleophile = "Nc1ccccc1"
    features = recommend.feat_molecular.featurize(electrophile, nucleophile)

    payload = recommend._compose_rule_features(
        "Ullmann_CN",
        features,
        role_pack=None,
        reactants=[electrophile, nucleophile],
    )

    assert payload["family"] == "Ullmann_CN"
    assert payload["electrophile"]["class"]
    assert payload["nucleophile"]["class"]
