"""Direction-neutral operator execution regressions."""

from __future__ import annotations

import pytest
from rdkit import Chem

from reactive_taxonomy import (
    BidirectionalReactionOperator,
    ReactionOperatorPrecedent,
    apply_forward_operator,
    apply_reaction_smarts,
    reverse_recovers_precursors,
)
from reactive_taxonomy.reaction_operators import operator_stereo_match


def test_matching_symmetric_spectators_preserves_opposite_assigned_stereo() -> None:
    spectators = ".".join(["CC(C)C"] * 18)
    assert operator_stereo_match("C[C@H](F)Br." + spectators, "C[C@@H](F)Br." + spectators) == "absent"


def test_partially_assigned_equivalent_centers_preserve_ambiguity_and_bound_search() -> None:
    a = "NC(F)C(Cl)[C@H](F)N"
    b = "NC(F)C(Cl)[C@@H](F)N"
    assert operator_stereo_match(a, b) == "connectivity_only"
    a = a.replace("(Cl)", "(C(F)C(Cl)C(F)C(Cl)C(F)C(Cl)C(F)O)")
    b = b.replace("(Cl)", "(C(F)C(Cl)C(F)C(Cl)C(F)C(Cl)C(F)O)")
    assert operator_stereo_match(a, b) == "unresolved"


def _operator() -> BidirectionalReactionOperator:
    precursor = "[#6:1]-[#6:2]-[#35].[#7&H3:3]"
    product = "[#6:1]-[#6:2]-[#7&H2:3]"
    return BidirectionalReactionOperator(
        operator_id="OP1:" + "a" * 64,
        realization_id="REAL2:" + "b" * 64,
        template_id="GRT3:test",
        abstraction_level="L2",
        forward_smarts=f"{precursor}>>{product}",
        reverse_smarts=f"{product}>>{precursor}",
        precursor_smarts=precursor,
        product_smarts=product,
        edit_tokens=(
            "formed:C-N:NONE>SINGLE",
            "hydrogen_change:N-H:SINGLE>NONE",
        ),
        operator_signature="test",
        stereo_policy="exact",
        observation_support=1,
        independent_reference_support=1,
        precedents=(
            ReactionOperatorPrecedent(
                reaction_id="reaction-1",
                reference_id="reference-1",
                precursor_smiles="CCBr.N",
                product_smiles="CCN",
            ),
        ),
    )


def test_forward_application_is_partner_order_invariant() -> None:
    operator = _operator()

    left = apply_forward_operator(operator, "CCBr.N")
    right = apply_forward_operator(operator, "N.CCBr")

    assert tuple(item.product_smiles for item in left) == ("CCN",)
    assert tuple(item.product_smiles for item in right) == ("CCN",)
    assert left[0].participating_precursor_smiles == "CCBr.N"
    assert right[0].participating_precursor_smiles == "CCBr.N"
    assert reverse_recovers_precursors(operator, left[0])
    assert left[0].mapped_product_smiles
    assert left[0].mapped_participating_precursor_smiles
    assert len(left[0].atom_correspondence) == 3
    assert {item.atom_map_number for item in left[0].atom_correspondence} <= {
        atom.GetAtomMapNum()
        for atom in Chem.MolFromSmiles(left[0].mapped_product_smiles).GetAtoms()
    }


def test_forward_application_ignores_unmatched_spectator_but_retains_indices() -> None:
    outcome = apply_forward_operator(_operator(), "CCBr.N.[Na+]")[0]

    assert outcome.product_smiles == "CCN"
    assert len(outcome.participating_component_indices) == 2
    assert "[Na+]" not in outcome.participating_precursor_smiles


def test_operator_serialization_is_deterministic() -> None:
    operator = _operator()

    restored = BidirectionalReactionOperator.from_dict(operator.to_dict())

    assert restored == operator
    assert restored.forward_operator_id == operator.forward_operator_id


def test_forward_application_can_use_virtual_copies_for_self_reaction() -> None:
    reaction_smarts = "[#6:1]-[#35].[#7&H2:2]>>[#6:1]-[#7&H1:2]"
    starting_material = "Brc1ccc(N)cc1"

    excluded = apply_reaction_smarts(reaction_smarts, starting_material)
    included = apply_reaction_smarts(
        reaction_smarts,
        starting_material,
        allow_self_reaction=True,
    )

    assert excluded == ()
    assert tuple(item.product_smiles for item in included) == (
        "Nc1ccc(Nc2ccc(Br)cc2)cc1",
    )
    outcome = included[0]
    assert outcome.assignment == (0, 0)
    assert outcome.participating_component_indices == (0,)
    assert outcome.reactant_stoichiometry == ((0, 2),)
    assert outcome.uses_virtual_copies is True
    assert outcome.participating_precursor_smiles == (
        "Nc1ccc(Br)cc1.Nc1ccc(Br)cc1"
    )
    assert {item.precursor_instance_index for item in outcome.atom_correspondence} == {
        0,
        1,
    }
    assert len(
        {item.atom_map_number for item in outcome.atom_correspondence}
    ) == len(outcome.atom_correspondence)
    assert "SELF_REACTION_USES_VIRTUAL_STOICHIOMETRIC_COPIES" in outcome.warnings

    explicit_copies = apply_reaction_smarts(
        reaction_smarts,
        f"{starting_material}.{starting_material}",
        allow_self_reaction=True,
    )
    assert explicit_copies
    assert all(not item.uses_virtual_copies for item in explicit_copies)


@pytest.mark.parametrize(("expected", "observed", "match"), [
    ("[13CH3]CO", "OCC", "absent"),
    ("[13CH3]CO", "C[13CH2]O", "absent"),
    ("[13CH3]CO", "OC[13CH3]", "exact"),
    ("[CH3:7][C@H:1](F)O", "C[C@H](F)O", "exact"),
    ("C[C@H](F)O", "C[C@@H](F)O", "absent"),
    ("C/C=C/C", "C/C=C\\C", "absent"),
    ("C/C=C/C", "CC=CC", "connectivity_only"),
    ("CC=CC", "C/C=C/C", "stereo_relaxed"),
    ("C[C@H](F)C[C@H](F)O", "C[C@@H](F)CC(F)O", "absent"),
    ("C[C@H](F)C[C@H](F)O", "C[C@H](F)CC(F)O", "connectivity_only"),
    ("CCO.[Na+].[Na+]", "CCO.[Na+]", "absent"),
    ("C(", "C(", "absent"),
])
def test_operator_stereo_match_preserves_constitution_and_requested_assignments(expected, observed, match):
    assert operator_stereo_match(expected, observed) == match
