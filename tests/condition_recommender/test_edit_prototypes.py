"""Tests for grammar-independent anonymous edit prototypes."""

from condition_recommender.edit_prototypes import (
    anonymous_edit_compatible,
    anonymous_edit_prototype,
    anonymous_edit_similarity,
)


def _signature(
    *,
    formed: list[str],
    broken: list[str],
    changed: list[str],
    ring_delta: int,
    ring_sizes: list[int],
) -> dict:
    return {
        "formed_bond_types": formed,
        "broken_bond_types": broken,
        "order_changes": changed,
        "event_count": 1,
        "topology": {
            "ring_count_delta": ring_delta,
            "formed_ring_sizes": ring_sizes,
        },
    }


def test_fischer_precursor_modes_share_anonymous_edit_neighborhood() -> None:
    separate_reactants = anonymous_edit_prototype(
        _signature(
            formed=["C-C:AROMATIC", "C-N:AROMATIC"],
            broken=["C-O:DOUBLE", "N-N:SINGLE"],
            changed=["C-C:SINGLE>AROMATIC"],
            ring_delta=1,
            ring_sizes=[5],
        )
    )
    preformed_hydrazone = anonymous_edit_prototype(
        _signature(
            formed=["C-C:AROMATIC", "C-N:AROMATIC"],
            broken=["C-N:DOUBLE", "N-N:SINGLE"],
            changed=["C-C:SINGLE>AROMATIC", "C-N:SINGLE>AROMATIC"],
            ring_delta=1,
            ring_sizes=[5],
        )
    )

    assert separate_reactants is not None
    assert preformed_hydrazone is not None
    assert anonymous_edit_compatible(
        separate_reactants,
        preformed_hydrazone,
    )
    assert (
        anonymous_edit_similarity(
            separate_reactants,
            preformed_hydrazone,
        )
        >= 0.6
    )


def test_unrelated_non_ring_edit_fails_prototype_compatibility_gate() -> None:
    annulation = anonymous_edit_prototype(
        _signature(
            formed=["C-C:SINGLE"],
            broken=["N-N:SINGLE"],
            changed=[],
            ring_delta=1,
            ring_sizes=[5],
        )
    )
    reduction = anonymous_edit_prototype(
        _signature(
            formed=[],
            broken=[],
            changed=["C-O:DOUBLE>SINGLE"],
            ring_delta=0,
            ring_sizes=[],
        )
    )

    assert annulation is not None
    assert reduction is not None
    assert not anonymous_edit_compatible(annulation, reduction)
    assert anonymous_edit_similarity(annulation, reduction) == 0.0
