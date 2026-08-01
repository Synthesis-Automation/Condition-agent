from types import SimpleNamespace

import pytest
from rdkit import Chem

from reactive_taxonomy.descriptors.activated_centers import (
    build_activated_center_context,
)
from reactive_taxonomy.descriptors.alkyl import build_alkyl_context
from reactive_taxonomy.descriptors.heteroatom import build_heteroatom_context
from reactive_taxonomy.descriptors.profiles import build_site_reactivity_profile
from reactive_taxonomy.descriptors.unsaturated import (
    build_alkenyl_context,
    build_alkynyl_context,
)


def _mol(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return molecule


def test_alkyl_substitution_and_beta_hydrogen_matrix() -> None:
    examples = (
        ("CBr", 0, "methyl", 0),
        ("CCBr", 1, "primary", 3),
        ("CC(Br)C", 1, "secondary", 6),
        ("CC(C)(C)Br", 1, "tertiary", 9),
        ("CC(C)(C)CBr", 4, "primary", 0),
    )
    for smiles, center, substitution, beta_hydrogens in examples:
        context = build_alkyl_context(_mol(smiles), center)[0]
        assert context.carbon_substitution == substitution
        assert context.beta_hydrogen_count == beta_hydrogens
    neopentyl = build_alkyl_context(_mol("CC(C)(C)CBr"), 4)[0]
    assert neopentyl.beta_branch_count == 2


def test_alkyl_cyclic_and_pi_adjacent_activation_are_orthogonal() -> None:
    cyclic = build_alkyl_context(_mol("BrC1CCCCC1"), 1)[0]
    benzylic = build_alkyl_context(_mol("BrCc1ccccc1"), 1)[0]
    allylic = build_alkyl_context(_mol("C=CCBr"), 2)[0]
    propargylic = build_alkyl_context(_mol("C#CCBr"), 2)[0]

    assert cyclic.cyclic and cyclic.ring_sizes == (6,)
    assert benzylic.benzylic and not benzylic.allylic
    assert allylic.allylic and not allylic.propargylic
    assert propargylic.propargylic and not propargylic.allylic


def test_alkenyl_and_alkynyl_contexts_cover_substitution_and_conjugation() -> None:
    ethene = build_alkenyl_context(_mol("C=C"), 0)[0]
    trans_butene = build_alkenyl_context(_mol("C/C=C/C"), 1)[0]
    cyclohexene = build_alkenyl_context(_mol("C1=CCCCC1"), 0)[0]
    diene = build_alkenyl_context(_mol("C=CC=C"), 0)[0]
    acetylene = build_alkynyl_context(_mol("C#C"), 0)[0]
    butyne = build_alkynyl_context(_mol("CC#CC"), 1)[0]
    phenylacetylene = build_alkynyl_context(
        _mol("c1ccccc1C#C"), 6
    )[0]

    assert ethene.alkene_class == "unsubstituted"
    assert trans_butene.endpoint_substitution == (1, 1)
    assert trans_butene.stereochemistry == "STEREOE"
    assert cyclohexene.cyclic and cyclohexene.ring_size == 6
    assert diene.conjugation_class == "extended_pi"
    assert acetylene.terminal
    assert not butyne.terminal
    assert phenylacetylene.conjugation_class == "aryl_conjugated"


def test_unsaturated_calculators_reject_nonmatching_centers() -> None:
    with pytest.raises(ValueError, match="double bond"):
        build_alkenyl_context(_mol("CCC"), 1)
    with pytest.raises(ValueError, match="triple bond"):
        build_alkynyl_context(_mol("C=C"), 0)


def test_activated_center_classification_matrix() -> None:
    examples = (
        ("CC(=O)Cl", 1, "acyl", "acyl_halide"),
        ("CC(=O)OC(=O)C", 1, "acyl", "anhydride"),
        ("CC(=O)OC", 1, "acyl", "ester_or_acid"),
        ("CC(=O)N", 1, "acyl", "amide"),
        ("CC=O", 1, "acyl", "aldehyde"),
        ("CC(=O)C", 1, "acyl", "ketone"),
        ("CS(=O)(=O)Cl", 1, "sulfonyl", "sulfonyl"),
        ("OP(=O)(O)O", 1, "phosphoryl", "phosphoryl"),
    )
    for smiles, center, kind, center_class in examples:
        context = build_activated_center_context(
            _mol(smiles), center
        )[0]
        assert context.context_kind == kind
        assert context.center_class == center_class

    acid_chloride = build_activated_center_context(
        _mol("CC(=O)Cl"), 1
    )[0]
    amide = build_activated_center_context(_mol("CC(=O)N"), 1)[0]
    assert acid_chloride.enolizable is True
    assert amide.enolizable is True


def test_heteroatom_context_separates_identity_resonance_and_group_burden() -> None:
    site = SimpleNamespace(details={})
    examples = (
        ("CN", 1, "primary", "amine_like", "localized"),
        ("CNC", 1, "secondary", "amine_like", "localized"),
        ("Nc1ccccc1", 0, "primary", "aryl_delocalized", "aryl_delocalized"),
        ("CC(=O)N", 3, "primary", "amide_like", "carbonyl_delocalized"),
        ("CS(=O)(=O)N", 4, "primary", "sulfonamide_like", "sulfonyl_delocalized"),
        ("CO", 1, "1_coordinate", "localized", "localized"),
        ("CS", 1, "1_coordinate", "localized", "localized"),
    )
    for smiles, center, substitution, lone_pair, resonance in examples:
        context = build_heteroatom_context(
            _mol(smiles), center, site
        )[0]
        assert context.substitution_class == substitution
        assert context.lone_pair_class == lone_pair
        assert context.resonance_class == resonance

    tert_butylamine = build_heteroatom_context(
        _mol("CC(C)(C)N"), 4, site
    )[0]
    assert tert_butylamine.attached_groups[0].attachment_carbon_class == (
        "tertiary"
    )
    assert tert_butylamine.alpha_branched_group_count == 1


def test_positive_charge_suppresses_heteroatom_lone_pair_availability() -> None:
    site = SimpleNamespace(details={})
    _, _, _, lone_pair_class, availability, acidity = (
        build_heteroatom_context(_mol("[NH3+]C"), 0, site)
    )
    assert lone_pair_class == "cationic"
    assert availability == "low"
    assert acidity == "not_applicable"


def test_negative_charge_is_retained_as_high_lone_pair_availability() -> None:
    site = SimpleNamespace(details={})
    _, _, _, lone_pair_class, availability, acidity = (
        build_heteroatom_context(_mol("[O-]C"), 0, site)
    )
    assert lone_pair_class == "anionic"
    assert availability == "high"
    assert acidity == "conjugate_base"


def test_unsupported_center_is_explicitly_unresolved() -> None:
    site = SimpleNamespace(
        hypothesis_id="manual:0",
        site_type="electrophilic_center",
        details={},
        confidence=1.0,
    )
    profile = build_site_reactivity_profile(
        _mol("[Na+]"),
        site,
        (),
        (),
        center_atom_index=0,
    )
    assert profile.context_kind == "other"
    assert profile.status == "unresolved"
    assert profile.context.reason == "unsupported_center_context"
