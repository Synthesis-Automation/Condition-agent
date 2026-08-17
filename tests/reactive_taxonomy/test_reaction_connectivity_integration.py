from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.reaction_edits import (
    normalize_inferred_scaffold_edits,
    normalize_mapped_edits,
    resolve_structural_evidence,
)
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


def _mapped(reaction: str):
    parsed = parse_reaction_smiles(reaction)
    assert parsed.valid
    result = normalize_mapped_edits(parsed.reactants, parsed.products)
    assert result.connectivity_edit_graph is not None
    return result


def test_mapped_multiple_bond_change_and_hydrogens_are_exact() -> None:
    result = _mapped("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    graph = result.connectivity_edit_graph
    assert graph is not None

    assert len(graph.bond_transitions) == 1
    transition = graph.bond_transitions[0]
    assert transition.before_state.order == "DOUBLE"
    assert transition.after_state.order == "SINGLE"
    assert transition.delta_units == -1
    assert transition.observation_scope == "observed_product"
    assert {
        (
            delta.atom.atom_map_number,
            delta.before_count,
            delta.after_count,
            delta.delta_count,
        )
        for delta in graph.hydrogen_deltas
    } == {(1, 2, 3, 1), (2, 2, 3, 1)}


def test_absent_leaving_group_is_projected_not_observed_no_bond() -> None:
    result = _mapped(
        "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]"
    )
    graph = result.connectivity_edit_graph
    assert graph is not None

    projected = next(
        transition
        for transition in graph.bond_transitions
        if {
            transition.atom_1.atom_map_number,
            transition.atom_2.atom_map_number,
        }
        == {1, 2}
    )
    formed = next(
        transition
        for transition in graph.bond_transitions
        if {
            transition.atom_1.atom_map_number,
            transition.atom_2.atom_map_number,
        }
        == {1, 3}
    )
    assert projected.before_state.order == "SINGLE"
    assert projected.after_state.state_kind == "endpoint_absent"
    assert projected.delta_units is None
    assert projected.observation_scope == "main_product_projection"
    assert formed.before_state.state_kind == "no_bond"
    assert formed.after_state.order == "SINGLE"
    assert formed.delta_units == 1
    assert formed.observation_scope == "observed_product"
    assert any(
        warning.startswith("PROJECTED_ATTACHMENT_NOT_FULLY_OBSERVED")
        for warning in graph.warnings
    )

    # The compatibility result remains unchanged during Phase 1.
    legacy = next(edit for edit in result.edits if edit.edit_type == "broken")
    assert legacy.evidence == "supplied_atom_mapping"
    assert legacy.old_order == "SINGLE"
    assert legacy.new_order is None


def test_unmapped_leaving_group_boundary_is_projected_when_product_is_complete() -> None:
    result = _mapped(
        "Br[CH2:1][CH3:2].[NH2:3][CH3:4]>>"
        "[CH3:2][CH2:1][NH:3][CH3:4]"
    )
    graph = result.connectivity_edit_graph
    assert graph is not None

    broken = next(edit for edit in result.edits if edit.edit_type == "broken")
    assert {broken.atom_1.element, broken.atom_2.element} == {"Br", "C"}
    assert broken.old_order == "SINGLE"
    projected = next(
        transition
        for transition in graph.bond_transitions
        if {transition.atom_1.element, transition.atom_2.element} == {"Br", "C"}
    )
    assert projected.after_state.state_kind == "endpoint_absent"
    assert projected.observation_scope == "main_product_projection"
    assert "PROJECTED_UNMAPPED_DEPARTING_BOUNDARY:1:Br" in result.warnings


def test_mapped_unmapped_leaving_group_and_unmapped_reaction_share_keys() -> None:
    mapped = featurize_reaction(
        "Br[CH2:1][CH3:2].[NH2:3][CH3:4]>>"
        "[CH3:2][CH2:1][NH:3][CH3:4]"
    )
    inferred = featurize_reaction("CCBr.NC>>CCNC")
    assert mapped.reaction_signature is not None
    assert inferred.reaction_signature is not None

    assert mapped.reaction_signature.formed_bond_types == ("C-N:SINGLE",)
    assert mapped.reaction_signature.broken_bond_types == ("Br-C:SINGLE",)
    assert (
        mapped.reaction_signature.transformation_signature_key
        == inferred.reaction_signature.transformation_signature_key
    )
    assert (
        mapped.reaction_signature.bond_edit_signature_key
        == inferred.reaction_signature.bond_edit_signature_key
    )
    assert (
        mapped.reaction_signature.handle_signature_key
        == inferred.reaction_signature.handle_signature_key
    )
    assert (
        mapped.reaction_signature.exact_signature_key
        == inferred.reaction_signature.exact_signature_key
    )


def test_mapped_mesylate_inventory_is_not_expanded_into_internal_edits() -> None:
    analysis = featurize_reaction(
        "[CH3:1][O:2][S:3]([CH3:4])(=[O:5])=[O:6]."
        "[NH2:7][CH3:8]>>[CH3:1][NH:7][CH3:8]"
    )

    signature = analysis.reaction_signature
    assert signature is not None
    assert signature.formed_bond_types == ("C-N:SINGLE",)
    assert signature.broken_bond_types == ("C-O:SINGLE",)
    assert signature.topology.ring_count_delta == 0
    assert (
        "REACTANT_ONLY_FRAGMENT_INTERNAL_BONDS_IGNORED:4"
        in analysis.warnings
    )


def test_mapped_tosylate_inventory_does_not_create_ring_deletion() -> None:
    analysis = featurize_reaction(
        "[CH3:1][O:2][S:3](=[O:4])(=[O:5])"
        "[c:6]1[cH:7][cH:8][cH:9][cH:10][cH:11]1."
        "[NH2:12][CH3:13]>>[CH3:1][NH:12][CH3:13]"
    )

    signature = analysis.reaction_signature
    assert signature is not None
    assert signature.formed_bond_types == ("C-N:SINGLE",)
    assert signature.broken_bond_types == ("C-O:SINGLE",)
    assert signature.topology.ring_count_delta == 0
    assert not any(
        token.startswith("C-C:") for token in signature.broken_bond_types
    )


def test_explicit_pair_addition_has_only_exact_bond_transitions() -> None:
    result = _mapped(
        "[CH2:1]=[CH2:2].[Br:3][Br:4]>>"
        "[CH2:1]([Br:3])[CH2:2][Br:4]"
    )
    graph = result.connectivity_edit_graph
    assert graph is not None

    assert len(graph.bond_transitions) == 4
    assert {
        transition.delta_units for transition in graph.bond_transitions
    } == {-1, 1}
    assert {
        transition.observation_scope for transition in graph.bond_transitions
    } == {"observed_product"}
    assert all(
        transition.after_state.state_kind != "endpoint_absent"
        for transition in graph.bond_transitions
    )


def test_formal_charge_change_is_retained_as_atom_state() -> None:
    result = _mapped("[NH3:1]>>[NH4+:1]")
    graph = result.connectivity_edit_graph
    assert graph is not None

    assert len(graph.atom_state_transitions) == 1
    transition = graph.atom_state_transitions[0]
    assert transition.reactant_atom.atom_map_number == 1
    assert transition.before_formal_charge == 0
    assert transition.after_formal_charge == 1
    assert transition.observation_scope == "observed_product"


def test_charge_only_mapping_survives_internal_normalization() -> None:
    parsed = parse_reaction_smiles("[Na:1]>>[Na+:1]")
    assert parsed.valid
    normalized = resolve_structural_evidence(
        parsed.reactants,
        parsed.products,
    )
    graph = normalized.connectivity_edit_graph
    assert graph is not None

    # Public edit resolution remains unchanged in Phase 1.
    assert not normalized.valid
    assert not normalized.edits
    assert len(graph.atom_state_transitions) == 1
    transition = graph.atom_state_transitions[0]
    assert transition.before_formal_charge == 0
    assert transition.after_formal_charge == 1


def test_product_only_atom_has_unknown_reactant_state() -> None:
    result = _mapped("[CH3:1]>>[CH3:1][Cl:2]")
    graph = result.connectivity_edit_graph
    assert graph is not None

    assert len(graph.bond_transitions) == 1
    transition = graph.bond_transitions[0]
    assert transition.before_state.state_kind == "unknown"
    assert transition.after_state.order == "SINGLE"
    assert transition.delta_units is None
    assert transition.observation_scope == "unresolved"
    assert any(
        warning.startswith("PRODUCT_ENDPOINT_WITHOUT_REACTANT_PROVENANCE")
        for warning in graph.warnings
    )


def test_partial_product_mapping_keeps_missing_endpoint_unknown() -> None:
    result = _mapped("[CH3:1][Br:2]>>C")
    graph = result.connectivity_edit_graph
    assert graph is not None

    assert len(graph.bond_transitions) == 1
    transition = graph.bond_transitions[0]
    assert transition.after_state.state_kind == "unknown"
    assert transition.observation_scope == "unresolved"
    assert transition.delta_units is None


def test_shadow_key_ignores_component_order_and_atom_map_values() -> None:
    forward = _mapped(
        "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]"
    )
    reordered = _mapped(
        "[NH2:30][CH3:40].[CH3:10][Br:20]>>"
        "[CH3:10][NH:30][CH3:40]"
    )
    assert forward.connectivity_edit_graph is not None
    assert reordered.connectivity_edit_graph is not None
    assert (
        forward.connectivity_edit_graph.shadow_key
        == reordered.connectivity_edit_graph.shadow_key
    )


def test_split_reagent_addition_uses_correspondence_scope() -> None:
    analysis = featurize_reaction("C=C.BrBr>>BrCCBr")
    assert analysis.observation is not None
    graph = analysis.observation.connectivity_edit_graph
    assert graph is not None
    assert analysis.evidence_quality == "global_atom_correspondence"
    assert {
        transition.observation_scope for transition in graph.bond_transitions
    } == {"correspondence_inference"}


def test_hydrosilane_correspondence_uses_hydrogen_delta_not_mapped_h() -> None:
    analysis = featurize_reaction(
        "C=C.C[SiH](C)C>>CC[Si](C)(C)C"
    )
    assert analysis.observation is not None
    graph = analysis.observation.connectivity_edit_graph
    assert graph is not None

    assert len(graph.hydrogen_deltas) == 2
    assert {delta.delta_count for delta in graph.hydrogen_deltas} == {-1, 1}
    assert all(
        delta.observation_scope == "correspondence_inference"
        for delta in graph.hydrogen_deltas
    )
    assert all(
        transition.atom_1.element != "H" and transition.atom_2.element != "H"
        for transition in graph.bond_transitions
    )


def test_inferred_correspondence_preserves_projected_attachment_scope() -> None:
    parsed = parse_reaction_smiles("CCBr.NC>>CCNC")
    assert parsed.valid
    normalized = normalize_inferred_scaffold_edits(
        parsed.reactants,
        parsed.products,
    )
    graph = normalized.connectivity_edit_graph
    assert graph is not None

    assert any(
        transition.after_state.state_kind == "endpoint_absent"
        and transition.observation_scope == "main_product_projection"
        for transition in graph.bond_transitions
    )
    assert any(
        transition.delta_units == 1
        and transition.observation_scope == "correspondence_inference"
        for transition in graph.bond_transitions
    )


def test_internal_dual_write_does_not_change_public_serialization() -> None:
    reaction = "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]"
    first = featurize_reaction(reaction)
    second = featurize_reaction(reaction)

    assert first.to_dict() == second.to_dict()
    assert first.reaction_signature is not None
    assert second.reaction_signature is not None
    assert (
        first.reaction_signature.signature_id
        == second.reaction_signature.signature_id
    )
    assert "connectivity_edit_graph" not in first.to_dict()
