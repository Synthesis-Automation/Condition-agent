"""Regression tests for conservative reaction-completeness assessment."""

from reactive_taxonomy import featurize_reaction


def test_structural_correspondence_verifies_product_atom_accounting() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )

    completeness = result.reaction_completeness
    assert completeness is not None
    assert completeness.status == "verified"
    assert completeness.product_heavy_atom_coverage == 1.0
    assert completeness.product_element_excess == {}
    assert completeness.reactant_element_excess == {
        "B": 1,
        "Br": 1,
        "O": 2,
    }
    assert completeness.warnings == ()
    assert result.reaction_signature is not None
    assert result.reaction_signature.completeness == completeness


def test_unaccounted_product_atoms_block_signature_generation() -> None:
    result = featurize_reaction(
        "[CH2:1]=[CH2:2]>>[CH3:1][CH2:2]C"
    )

    completeness = result.reaction_completeness
    assert completeness is not None
    assert completeness.status == "incomplete"
    assert completeness.product_element_excess == {"C": 1}
    assert completeness.product_heavy_atom_coverage == 0.666667
    assert completeness.suspected_missing_reactant
    assert "UNACCOUNTED_PRODUCT_HEAVY_ATOMS" in completeness.warnings
    assert result.reaction_signature is None


def test_exact_insufficient_partner_multiplicity_is_explicitly_inferred() -> None:
    result = featurize_reaction(
        "Sc1ccccc1.O=[N+]([O-])c1c(F)cccc1F"
        ">>O=[N+]([O-])c1c(Sc2ccccc2)cccc1Sc1ccccc1",
        label_style="ascii",
    )

    completeness = result.reaction_completeness
    assert completeness is not None
    assert completeness.status == "verified"
    assert completeness.suspected_insufficient_reactant_multiplicity
    assert not completeness.suspected_missing_reactant
    assert "INFERRED_REACTANT_MULTIPLICITY" in completeness.warnings
    assert result.reaction_signature is not None
    assert result.reaction_signature.event_count == 2
    assert len(result.reactants) == 3
    assert result.reactants[2].inferred_copy_of_component_index == 0


def test_partial_mapping_is_reported_per_heavy_atom() -> None:
    result = featurize_reaction("[CH3:1][OH]>>[CH2:1]=O")

    completeness = result.reaction_completeness
    assert completeness is not None
    assert completeness.status == "unresolved"
    assert completeness.reactant_mapping_coverage == 0.5
    assert completeness.product_mapping_coverage == 0.5
    assert "PARTIAL_ATOM_MAPPING" in completeness.warnings
    assert "PARTIAL_ATOM_MAPPING" in result.warnings


def test_product_map_numbers_absent_from_reactants_are_reported() -> None:
    result = featurize_reaction("[CH3:1][OH:2]>>[CH2:1]=[O:3]")

    completeness = result.reaction_completeness
    assert completeness is not None
    assert completeness.status == "unresolved"
    assert completeness.shared_mapped_heavy_atom_count == 1
    assert "PRODUCT_MAPS_MISSING_FROM_REACTANTS" in completeness.warnings
    assert "REACTANT_MAPS_MISSING_FROM_PRODUCTS" in completeness.warnings


def test_unresolved_atom_provenance_remains_explicit() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1"
    )

    completeness = result.reaction_completeness
    assert completeness is not None
    assert completeness.status == "unresolved"
    assert completeness.product_heavy_atom_coverage is None
    assert completeness.warnings == ("REACTION_COMPLETENESS_UNRESOLVED",)
