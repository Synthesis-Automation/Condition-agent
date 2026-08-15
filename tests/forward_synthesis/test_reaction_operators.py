"""Direction-neutral operator execution regressions."""

from __future__ import annotations

from rdkit import Chem

from reactive_taxonomy import (
    BidirectionalReactionOperator,
    ReactionOperatorPrecedent,
    apply_forward_operator,
    reverse_recovers_precursors,
)


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
