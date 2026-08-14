"""Tests for interpretation-independent anonymous edit prototypes."""

from dataclasses import asdict

from reactive_taxonomy import featurize_reaction

from condition_recommender.edit_prototypes import (
    anonymous_edit_compatible,
    anonymous_edit_prototype,
    anonymous_edit_prototype_from_hypothesis,
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


def test_ambiguous_fischer_hypotheses_share_one_anonymous_prototype() -> None:
    analysis = featurize_reaction(
        "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
        ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
    )

    prototypes = tuple(
        anonymous_edit_prototype_from_hypothesis(asdict(hypothesis))
        for hypothesis in analysis.edit_hypotheses
    )

    assert len(prototypes) == 2
    assert all(prototype is not None for prototype in prototypes)
    assert prototypes[0] == prototypes[1]
    assert prototypes[0].formed_element_pairs == ("C-C", "C-N")
    assert prototypes[0].broken_element_pairs == ("C-O", "N-N")


def test_center_state_separates_alkylation_from_amide_formation() -> None:
    mapped_mesylate = featurize_reaction(
        "[CH3:1][O:2][S:3]([CH3:4])(=[O:5])=[O:6]."
        "[NH2:7][CH3:8]>>[CH3:1][NH:7][CH3:8]"
    )
    partially_mapped_mesylate = featurize_reaction(
        "CS(=O)(=O)O[CH3:1].[NH2:7][CH3:8]"
        ">>[CH3:1][NH:7][CH3:8]"
    )
    amidation = featurize_reaction(
        "[CH3:1][C:2](=[O:3])[O:4].[NH2:5][CH3:6]"
        ">>[CH3:1][C:2](=[O:3])[NH:5][CH3:6]"
    )
    assert mapped_mesylate.reaction_signature is not None
    assert partially_mapped_mesylate.reaction_signature is not None
    assert amidation.reaction_signature is not None

    mapped_prototype = anonymous_edit_prototype(
        asdict(mapped_mesylate.reaction_signature)
    )
    partial_prototype = anonymous_edit_prototype(
        asdict(partially_mapped_mesylate.reaction_signature)
    )
    amide_prototype = anonymous_edit_prototype(
        asdict(amidation.reaction_signature)
    )
    assert mapped_prototype is not None
    assert partial_prototype is not None
    assert amide_prototype is not None

    assert mapped_prototype.formed_atom_state_pairs == ("C@SP3-N",)
    assert amide_prototype.formed_atom_state_pairs == ("C@SP2-N",)
    assert anonymous_edit_compatible(mapped_prototype, partial_prototype)
    assert anonymous_edit_similarity(mapped_prototype, partial_prototype) >= 0.6
    assert not anonymous_edit_compatible(mapped_prototype, amide_prototype)
    assert anonymous_edit_similarity(mapped_prototype, amide_prototype) == 0.0
