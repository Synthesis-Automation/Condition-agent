from pathlib import Path

from condition_recommender import recommend_generic_conditions
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.generic_retrieval import (
    generic_signature_similarity,
    load_generic_retrieval_rules,
    retrieve_compatible_generic_pool,
    retrieve_generic_pool,
)


def _signature(
    token: str,
    *,
    exact: str = "exact-a",
    handles: str = "handles-a",
    transformation: str = "transformation-a",
    bond: str = "bond-a",
    environment: str = "environment-a",
    family: str | None = "family-a",
    family_confidence: float = 1.0,
) -> dict:
    return {
        "signature_id": token,
        "exact_signature_key": exact,
        "handle_signature_key": handles,
        "transformation_signature_key": transformation,
        "bond_edit_signature_key": bond,
        "environment_signature_key": environment,
        "formed_bond_types": ["C-N:SINGLE"],
        "broken_bond_types": ["B-C:SINGLE"],
        "order_changes": [],
        "partners": [
            {
                "handle_tokens": ["B(OH)2"],
                "anchor_contexts": ["Ar"],
                "steric": {"class": "open"},
                "electronic": {"class": "neutral"},
                "nearby_groups": [],
            },
            {
                "handle_tokens": ["N-H"],
                "anchor_contexts": ["Ar"],
                "steric": {"class": "open"},
                "electronic": {"class": "neutral"},
                "nearby_groups": [],
            },
        ],
        "spectator_groups": [],
        "transformation_class": "generic_c_n_coupling",
        "named_family": family,
        "family_confidence": family_confidence,
    }


def _record(index: int, signature: dict, *, tier: str = "verified") -> dict:
    return {
        "admission_tier": tier,
        "reaction_id": f"reaction-{index}",
        "observation_id": f"observation-{index}",
        "reaction_smiles": "C.N>>CN",
        "yield_pct": 60.0 + index,
        "source_dataset": f"dataset-{index % 2}",
        "reaction_signature": signature,
        "resolved_recipe_id": f"RCR1:{index % 2}",
        "resolved_recipe": {"recipe_id": f"RCR1:{index % 2}"},
        "condition_resolution": {"has_uncertainty": False},
    }


def test_generic_index_admits_only_usable_verified_records() -> None:
    index = build_generic_index(
        [
            _record(1, _signature("one")),
            _record(2, _signature("two"), tier="review"),
            {**_record(3, _signature("three")), "resolved_recipe": None},
        ]
    )
    assert len(index.rows) == 1
    assert index.exact["exact-a"] == (0,)


def test_retrieval_uses_exact_signature_when_supported() -> None:
    query = _signature("query")
    records = [_record(index, _signature(str(index))) for index in range(3)]
    records.append(
        _record(9, _signature("other", exact="exact-b", handles="handles-b"))
    )
    level, pool = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "exact_signature"
    assert len(pool) == 3


def test_retrieval_falls_back_exact_to_handles_then_family() -> None:
    query = _signature("query")
    records = [
        _record(
            index,
            _signature(str(index), exact=f"other-{index}"),
        )
        for index in range(3)
    ]
    level, _ = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "handle_signature"

    family_records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"other-{index}",
                handles=f"other-{index}",
                transformation=f"other-{index}",
            ),
        )
        for index in range(3)
    ]
    level, pool = retrieve_generic_pool(query, build_generic_index(family_records))
    assert level == "named_family"
    assert len(pool) == 3


def test_low_confidence_family_uses_transformation_fallback() -> None:
    query = _signature("query", family_confidence=0.5)
    records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"exact-{index}",
                handles=f"handles-{index}",
                family="different-family",
            ),
        )
        for index in range(3)
    ]
    level, pool = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "transformation_signature"
    assert len(pool) == 3


def test_family_fallback_never_crosses_bond_edit_gate() -> None:
    query = _signature("query")
    records = [
        _record(
            index,
            _signature(
                str(index),
                exact=f"exact-{index}",
                handles=f"handles-{index}",
                transformation=f"transformation-{index}",
                bond="incompatible-bond",
            ),
        )
        for index in range(3)
    ]
    level, pool = retrieve_generic_pool(query, build_generic_index(records))
    assert level == "no_compatible_bond_edit"
    assert pool == ()


def test_compatibility_exclusion_continues_to_relaxed_tier() -> None:
    query = _signature("query")
    query["spectator_groups"] = [
        {"group_id": "aldehyde", "tags": ["oxidation_sensitive"]}
    ]
    unsafe = _record(0, _signature("unsafe"))
    unsafe["resolved_recipe"] = {
        "recipe_id": unsafe["resolved_recipe_id"],
        "oxidants": [
            {
                "roles": [
                    {"family_id": "peroxides_hydroperoxides", "role_id": "oxidant"}
                ]
            }
        ],
    }
    safe_records = [
        _record(
            index,
            _signature(str(index), exact=f"safe-exact-{index}"),
        )
        for index in range(1, 4)
    ]
    level, accepted, raw_count, excluded_count = retrieve_compatible_generic_pool(
        query,
        build_generic_index([unsafe, *safe_records]),
        minimum_pool_size=2,
    )
    assert level == "handle_signature"
    assert len(accepted) == 3
    assert raw_count == 4
    assert excluded_count == 1


def test_similarity_does_not_treat_missing_features_as_matches() -> None:
    score, components = generic_signature_similarity(
        {
            "formed_bond_types": ["C-N:SINGLE"],
            "partners": [],
        },
        {
            "formed_bond_types": ["C-N:SINGLE"],
            "partners": [],
        },
    )
    assert components["edit_topology"] == 1.0
    assert components["handles"] == 0.0
    assert components["environment"] == 0.0
    assert 0.0 < score < 1.0


def test_generic_retrieval_weights_are_normalized() -> None:
    rules = load_generic_retrieval_rules()
    assert round(sum(rules["similarity_weights"].values()), 10) == 1.0
    assert round(sum(rules["ranking_weights"].values()), 10) == 1.0


def test_real_pilot_returns_resolved_recipe() -> None:
    path = Path("results/generic_conversion_chan_lam_pilot/records.jsonl")
    result = recommend_generic_conditions(
        "c1ccc2[nH]cnc2c1.COc1ccc(B(O)O)cc1>>"
        "COc1ccc(-n2cnc3ccccc32)cc1",
        records_path=path,
        top_k=3,
        minimum_pool_size=1,
    )
    assert result.valid
    assert result.retrieval_level == "exact_signature"
    assert result.candidate_count == 1
    assert result.compatible_candidate_count == 1
    assert result.excluded_candidate_count == 0
    assert result.recommendations
    assert result.recommendations[0].recipe_id.startswith("RCR1:")
    assert result.recommendations[0].resolved_recipe["recipe_id"].startswith("RCR1:")
    assert 0.0 <= result.recommendations[0].compatibility_score <= 1.0
