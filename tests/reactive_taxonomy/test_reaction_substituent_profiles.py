"""Graph-derived R-group profile and core-event relationship regressions."""

from __future__ import annotations

from reactive_taxonomy import (
    build_substituent_profile,
    featurize_reaction,
)
from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles
from reactive_taxonomy.reaction_display_projection import (
    build_reaction_display_projection,
)
from reactive_taxonomy.reaction_render_context import (
    reaction_render_context_from_analysis,
)


def _oxygen_bound_profile(smiles: str):
    molecule = parse_smiles(smiles)
    assert molecule is not None
    oxygen = next(
        atom for atom in molecule.GetAtoms() if atom.GetSymbol() == "O"
    )
    attachment = next(
        atom for atom in oxygen.GetNeighbors() if atom.GetAtomicNum() > 1
    )
    fragment = tuple(
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if atom.GetIdx() != oxygen.GetIdx()
    )
    bond = molecule.GetBondBetweenAtoms(oxygen.GetIdx(), attachment.GetIdx())
    return build_substituent_profile(
        molecule,
        fragment_atom_indices=fragment,
        attachment_atom_index=attachment.GetIdx(),
        core_atom_index=oxygen.GetIdx(),
        attachment_bond_order=str(bond.GetBondType()).upper(),
    )


def test_alcohol_attachment_profiles_distinguish_substitution_degree() -> None:
    methyl = _oxygen_bound_profile("CO")
    primary = _oxygen_bound_profile("CCO")
    secondary = _oxygen_bound_profile("CC(C)O")
    tertiary = _oxygen_bound_profile("CC(C)(C)O")

    assert methyl.carbon_substitution == "methyl"
    assert primary.carbon_substitution == "primary"
    assert secondary.carbon_substitution == "secondary"
    assert tertiary.carbon_substitution == "tertiary"
    assert {value.base_class for value in (methyl, primary, secondary, tertiary)} == {
        "alkyl"
    }


def test_attachment_profiles_capture_resonance_and_ring_context() -> None:
    phenyl = _oxygen_bound_profile("Oc1ccccc1")
    benzyl = _oxygen_bound_profile("OCC1=CC=CC=C1")
    allyl = _oxygen_bound_profile("OCC=C")
    propargyl = _oxygen_bound_profile("OCC#C")
    cyclohexyl = _oxygen_bound_profile("OC1CCCCC1")

    assert phenyl.base_class == "aryl"
    assert phenyl.carbon_substitution == "not_applicable"
    assert benzyl.benzylic and benzyl.carbon_substitution == "primary"
    assert allyl.allylic and not allyl.propargylic
    assert propargyl.propargylic and not propargyl.allylic
    assert cyclohexyl.base_class == "ring_aliphatic"
    assert cyclohexyl.cyclic and cyclohexyl.ring_sizes == (6,)


def test_profile_identity_is_invariant_to_equivalent_smiles_order() -> None:
    forward = _oxygen_bound_profile("CCCO")
    reversed_smiles = _oxygen_bound_profile("OCCC")

    assert forward.profile_id == reversed_smiles.profile_id
    assert forward.feature_tokens == reversed_smiles.feature_tokens


def test_aryl_profiles_preserve_ortho_meta_para_substituents() -> None:
    reactions = {
        "ortho": "Brc1ccccc1F.CN>>CNc1ccccc1F",
        "meta": "Brc1cccc(F)c1.CN>>CNc1cccc(F)c1",
        "para": "Brc1ccc(F)cc1.CN>>CNc1ccc(F)cc1",
    }
    observed = {}
    for expected, reaction in reactions.items():
        analysis = featurize_reaction(reaction)
        assert analysis.reaction_core is not None
        relations = tuple(
            relation
            for remote in analysis.reaction_core.remote_subgraphs
            if remote.side == "reactant"
            for port in remote.attachment_ports
            for relation in port.substituent_profile.aromatic_substituent_relations
            if relation.substituent_fragment_smiles == "F"
        )
        assert relations
        assert {value.positional_relation for value in relations} == {expected}
        observed[expected] = {
            port.substituent_profile.profile_id
            for remote in analysis.reaction_core.remote_subgraphs
            if remote.side == "reactant"
            for port in remote.attachment_ports
            if port.substituent_profile.aromatic_substituent_relations
        }

    assert len({next(iter(values)) for values in observed.values()}) == 3


def test_core_and_display_serialize_the_same_profile_contract() -> None:
    analysis = featurize_reaction("CC(=O)O.CN>>CC(=O)NC")
    assert analysis.reaction_core is not None
    profiles = tuple(
        port.substituent_profile
        for remote in analysis.reaction_core.remote_subgraphs
        for port in remote.attachment_ports
    )

    assert profiles
    assert all(profile.feature_tokens for profile in profiles)
    assert all(profile.definition_version == "substituent_profiles.v1" for profile in profiles)
    display = build_reaction_display_projection(
        reaction_render_context_from_analysis(analysis)
    )
    assert display.substituents
    assert all(
        value.substituent_profile.definition_version == "substituent_profiles.v1"
        for value in display.substituents
    )


def test_multisite_core_records_same_molecule_shortest_path() -> None:
    analysis = featurize_reaction(
        "Cc1c([N+](=O)[O-])cnc2c1c("
        "C1=CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
        ">>Cc1c(N)cnc2c1c("
        "C1CCN(C(=O)OC(C)(C)C)C(C)(C)C1)cn2C"
    )
    assert analysis.reaction_core is not None
    assert analysis.reaction_core.event_count == 2
    assert len(analysis.reaction_core.event_relations) == 1
    relation = analysis.reaction_core.event_relations[0]

    assert relation.relation_type == "same_component"
    assert relation.shared_reactant_component_indices == (0,)
    assert any(
        path.side == "reactant" and path.bond_count > 0
        for path in relation.shortest_paths
    )


def test_patterns_and_text_cite_core_events_and_r_profiles() -> None:
    analysis = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")
    assert analysis.interpretation is not None
    match = next(
        value
        for value in analysis.interpretation.pattern_matches
        if value.pattern_id == "sp2_c_n_substitution_like"
    )
    assert match.matched_core_event_ids
    assert match.matched_substituent_profile_ids
    assert match.covered_core_event_fraction == 1.0
    assert analysis.reaction_label is not None
    assert analysis.reaction_label.core_event_ids == match.matched_core_event_ids
    assert set(match.matched_substituent_profile_ids).issubset(
        analysis.reaction_label.substituent_profile_ids
    )


def test_core_event_and_profile_provenance_is_partner_order_invariant() -> None:
    forward = featurize_reaction("Brc1ccccc1.CN>>CNc1ccccc1")
    reversed_order = featurize_reaction("CN.Brc1ccccc1>>CNc1ccccc1")
    assert forward.reaction_core is not None
    assert reversed_order.reaction_core is not None

    assert tuple(event.event_id for event in forward.reaction_core.events) == tuple(
        event.event_id for event in reversed_order.reaction_core.events
    )
    assert {
        port.substituent_profile.profile_id
        for remote in forward.reaction_core.remote_subgraphs
        for port in remote.attachment_ports
    } == {
        port.substituent_profile.profile_id
        for remote in reversed_order.reaction_core.remote_subgraphs
        for port in remote.attachment_ports
    }
    assert forward.reaction_label is not None
    assert reversed_order.reaction_label is not None
    assert forward.reaction_label.core_event_ids == (
        reversed_order.reaction_label.core_event_ids
    )
