from condition_recommender.signature_features import (
    environment_profile_similarity,
)


def _profile(
    *,
    kind: str = "aromatic",
    access: str = "open",
    steric_score: float = 0.1,
    axis: str = "electronic_demand",
    electronic_class: str = "balanced",
    electronic_score: float = 0.0,
) -> dict:
    return {
        "context_kind": kind,
        "context": {
            "context_kind": kind,
            "ring_family": "benzene",
            "ring_sizes": [6],
            "ortho_occupancy_count": 0,
            "ortho_capacity": 2,
            "ortho_burden_class": "none",
            "heteroatoms": [],
        },
        "steric": {
            "accessibility_class": access,
            "accessibility_score": steric_score,
        },
        "electronic": {
            "activation_axis": axis,
            "activation_class": electronic_class,
            "activation_score": electronic_score,
        },
        "reactive_center": {},
        "modifiers": [],
    }


def _signature(profile: dict, role: str = "electrophile") -> dict:
    return {
        "partners": [
            {
                "role": role,
                "reactivity_profile": profile,
                "nearby_groups": [],
            }
        ]
    }


def test_profile_similarity_uses_numeric_distance_within_the_same_axis() -> None:
    query = _signature(_profile(electronic_score=0.1))
    close = _signature(_profile(electronic_score=0.15))
    distant = _signature(
        _profile(
            electronic_class="electron_poor",
            electronic_score=0.9,
        )
    )
    assert environment_profile_similarity(
        query, close
    ) > environment_profile_similarity(query, distant)


def test_profile_similarity_does_not_compare_scores_across_axes() -> None:
    query = _signature(_profile(electronic_score=0.1))
    same_axis = _signature(_profile(electronic_score=0.1))
    other_axis = _signature(
        _profile(
            axis="lone_pair_availability",
            electronic_class="balanced",
            electronic_score=0.1,
        )
    )
    assert environment_profile_similarity(
        query, same_axis
    ) > environment_profile_similarity(query, other_axis)


def test_profile_similarity_aligns_partners_by_role() -> None:
    query = _signature(_profile(), role="electrophile")
    wrong_role = _signature(_profile(), role="nucleophile")
    assert environment_profile_similarity(query, wrong_role) == 0.0
