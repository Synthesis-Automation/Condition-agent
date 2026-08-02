"""Regressions for incomplete but structurally informative reactions."""

from reactive_taxonomy import featurize_reaction


def test_missing_chlorine_source_yields_partial_acyl_substitution() -> None:
    result = featurize_reaction(
        "O=C(O)c1cccc(I)c1>>O=C(Cl)c1cccc(I)c1"
    )

    assert result.valid
    assert result.reaction_label.text == "R–C(=O)–O → R–C(=O)–Cl"
    assert result.reaction_label.status == "partial_product_correspondence"
    assert result.evidence_quality == "partial_product_correspondence"
    assert result.transformation_class == "acyl_heteroatom_substitution"
    assert result.named_family is None
    assert result.reaction_signature is None
    assert result.reaction_completeness.status == "incomplete"
    assert result.reaction_label.status == "partial_product_correspondence"
    assert result.reaction_label.basis == "partial_product_correspondence"
    assert "Ar–I" not in result.reaction_label.text

    observation = result.partial_product_transformation
    assert observation is not None
    assert observation.transformation_type == "attachment_replacement"
    assert observation.reactant_center.element == "C"
    assert observation.removed_attachment.element == "O"
    assert observation.added_attachment.element == "Cl"
    assert observation.missing_product_atom_elements == ("Cl",)
    assert observation.product_heavy_atom_coverage == 0.9
    assert "PRODUCT_ATOM_SOURCE_UNRESOLVED:Cl" in result.warnings
    assert "PARTIAL_PRODUCT_CORRESPONDENCE" in result.warnings
    descriptor = result.fallback_descriptor
    assert descriptor is not None
    assert descriptor.partial_transformation_key.startswith("PTS1:")
    assert len(descriptor.source_requirements) == 1
    requirement = descriptor.source_requirements[0]
    assert requirement.requirement_id.startswith("FSR1:")
    assert requirement.element_counts == {"Cl": 1}
    assert requirement.attachment_element == "Cl"
    assert requirement.removed_attachment_element == "O"


def test_partial_attachment_replacement_is_not_acyl_specific() -> None:
    result = featurize_reaction("CC(C)O>>CC(C)Cl")

    assert result.reaction_label.text == "C–O → C–Cl"
    assert result.transformation_class == "attachment_replacement"
    assert result.partial_product_transformation is not None
    assert result.reaction_signature is None
    assert result.reaction_label.basis == "partial_product_correspondence"


def test_partial_replacement_ascii_label_is_deterministic() -> None:
    result = featurize_reaction(
        "CC(=O)O>>CC(=O)F",
        label_style="ascii",
    )

    assert result.reaction_label.text == "R-C(=O)-O -> R-C(=O)-F"
    assert result.partial_product_transformation is not None


def test_product_growth_without_attachment_exchange_is_not_forced() -> None:
    result = featurize_reaction("CC>>CCC")

    assert result.partial_product_transformation is None
    assert result.reaction_signature is None


def test_missing_cyanide_source_retains_rooted_product_fragment() -> None:
    result = featurize_reaction("Brc1ccccc1>>N#Cc1ccccc1")

    assert result.reaction_label.text == "C–Br → C–C≡N"
    assert result.reaction_label.status == "partial_product_correspondence"
    assert result.reaction_signature is None

    observation = result.partial_product_transformation
    assert observation is not None
    assert observation.transformation_type == "attachment_fragment_replacement"
    assert observation.removed_fragment_smiles == "Br"
    assert observation.missing_product_atom_elements == ("C", "N")
    assert observation.installed_fragment.canonical_fragment_smiles == "C#N"
    assert observation.installed_fragment.rooted_fragment_smiles == "*C#N"
    assert observation.installed_fragment.internal_bond_types == (
        "C-N:TRIPLE",
    )
    assert observation.installed_fragment.element_counts == {"C": 1, "N": 1}
    assert observation.installed_fragment.source_status == "unresolved"
    assert observation.installed_fragment.source_candidates == ()
    assert observation.installed_fragment.atom_references[0].side == "product"
    assert len(observation.product_atom_provenance) == 8
    assert sum(
        item.source_kind == "reactant_correspondence"
        for item in observation.product_atom_provenance
    ) == 6
    unresolved = tuple(
        item
        for item in observation.product_atom_provenance
        if item.source_kind == "unresolved"
    )
    assert len(unresolved) == 2
    assert all(item.source_atom is None for item in unresolved)
    assert "PRODUCT_FRAGMENT_SOURCE_UNRESOLVED" in result.warnings


def test_missing_azide_source_retains_multi_atom_origin_gap() -> None:
    result = featurize_reaction(
        "Brc1ccccc1>>[N-]=[N+]=Nc1ccccc1",
        label_style="ascii",
    )

    observation = result.partial_product_transformation
    assert observation is not None
    assert result.reaction_label.text == "C-Br -> C-[N-]=[N+]=N"
    assert observation.installed_fragment.element_counts == {"N": 3}
    assert observation.installed_fragment.internal_bond_types == (
        "N-N:DOUBLE",
        "N-N:DOUBLE",
    )
    assert observation.installed_fragment.rooted_fragment_smiles == (
        "*N=[N+]=[N-]"
    )
    assert len(observation.installed_fragment.atom_references) == 3
    assert len(observation.installed_fragment.attachments) == 1


def test_product_fragment_key_is_serialization_invariant() -> None:
    left_first = featurize_reaction(
        "Brc1ccccc1>>[N-]=[N+]=Nc1ccccc1"
    )
    ring_first = featurize_reaction(
        "Brc1ccccc1>>c1ccccc1N=[N+]=[N-]"
    )

    left_observation = left_first.partial_product_transformation
    ring_observation = ring_first.partial_product_transformation
    assert left_observation is not None
    assert ring_observation is not None
    assert (
        left_observation.installed_fragment.fragment_key
        == ring_observation.installed_fragment.fragment_key
    )


def test_middle_side_cyanide_is_structural_support_not_invented_mapping() -> None:
    result = featurize_reaction(
        "Brc1ccccc1>[Na+].[C-]#N>N#Cc1ccccc1",
        label_style="ascii",
    )

    observation = result.partial_product_transformation
    assert observation is not None
    fragment = observation.installed_fragment
    assert fragment.source_status == "agent_supported"
    assert len(fragment.source_candidates) == 1
    assert fragment.source_candidates[0].side == "agent"
    assert fragment.source_candidates[0].canonical_smiles == "[C-]#N"
    assert fragment.source_candidates[0].evidence == (
        "agent_fragment_subgraph_support"
    )
    assert all(
        reference.side == "product"
        for reference in fragment.atom_references
    )
    agent_supported = tuple(
        item
        for item in observation.product_atom_provenance
        if item.source_kind == "agent_supported"
    )
    assert len(agent_supported) == 2
    assert all(item.source_atom is None for item in agent_supported)
    assert result.reaction_label.text == "C-Br -> C-C#N"
    assert "PRODUCT_FRAGMENT_SOURCE_AGENT_SUPPORTED" in result.warnings
    assert "PRODUCT_FRAGMENT_SOURCE_UNRESOLVED" not in result.warnings


def test_multiple_structural_agent_sources_preserve_ambiguity() -> None:
    result = featurize_reaction(
        "Brc1ccccc1>[C-]#N.C[Si](C)(C)C#N>N#Cc1ccccc1"
    )

    observation = result.partial_product_transformation
    assert observation is not None
    fragment = observation.installed_fragment
    assert fragment.source_status == "ambiguous"
    assert len(fragment.source_candidates) == 2
    assert "AMBIGUOUS_PRODUCT_FRAGMENT_SOURCES" in result.warnings
    assert result.reaction_label.text == "C–Br → C–C≡N"


def test_partial_nitrile_label_preserves_ascii_triple_bond_style() -> None:
    result = featurize_reaction(
        "Brc1ccccc1>>N#Cc1ccccc1",
        label_style="ascii",
    )

    assert result.reaction_label.text == "C-Br -> C-C#N"


def test_distinct_installed_fragments_receive_distinct_graph_keys() -> None:
    cyanide = featurize_reaction("Brc1ccccc1>>N#Cc1ccccc1")
    ethyl = featurize_reaction("Brc1ccccc1>>CCc1ccccc1")

    cyanide_observation = cyanide.partial_product_transformation
    ethyl_observation = ethyl.partial_product_transformation
    assert cyanide_observation is not None
    assert ethyl_observation is not None
    assert (
        cyanide_observation.installed_fragment.fragment_key
        != ethyl_observation.installed_fragment.fragment_key
    )


def test_general_fragment_replacement_supports_multi_atom_leaving_branch() -> None:
    result = featurize_reaction(
        "O=S(=O)(Oc1ccccc1)C(F)(F)F>>N#Cc1ccccc1",
        label_style="ascii",
    )

    observation = result.partial_product_transformation
    assert observation is not None
    assert len(observation.removed_fragment_atom_indices) > 1
    assert observation.removed_fragment_smiles == "O=S(=O)(O)C(F)(F)F"
    assert observation.installed_fragment.rooted_fragment_smiles == "*C#N"
