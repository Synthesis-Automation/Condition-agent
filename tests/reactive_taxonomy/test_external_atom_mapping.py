"""Tests for optional external atom-mapping validation and provenance."""

from __future__ import annotations

import builtins
from dataclasses import replace

from reactive_taxonomy import (
    AtomMappingProviderMetadata,
    EXTERNAL_MAPPING_EVIDENCE,
    RxnMapperProvider,
    analyze_reaction_with_external_mapping,
    build_reaction_display_projection,
    build_reaction_review_summary,
    featurize_reaction,
    reaction_render_context_from_analysis,
    validate_external_atom_mapping,
)
from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


FISCHER_REACTION = (
    "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
    ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
)
FISCHER_RXNMAPPER_OUTPUT = (
    "Cl.N[NH:6][c:5]1[cH:4][cH:3][c:2]([F:1])[cH:10][cH:9]1."
    "O=[C:7]1[CH2:8][CH2:11][CH2:12][CH2:13][CH2:14]1"
    ">>[F:1][c:2]1[cH:3][cH:4][c:5]2[nH:6][c:7]3"
    "[c:8]([c:9]2[cH:10]1)[CH2:11][CH2:12][CH2:13][CH2:14]3"
)
METADATA = AtomMappingProviderMetadata(
    provider_id="rxnmapper",
    provider_version="test",
    model_id="fixture",
    model_sha256="abc",
)

RESOLVED_SUZUKI_REACTION = (
    "CC(C)(C)c1ccc(OS(C)(=O)=O)cc1.OB(O)c1ccccc1"
    ">>CC(C)(C)c1ccc(-c2ccccc2)cc1"
)
RESOLVED_SUZUKI_MAPPING = (
    "CS(=O)(=O)[O:17][c:8]1[cH:7][cH:6][c:5]"
    "([C:2]([CH3:1])([CH3:3])[CH3:4])[cH:16][cH:15]1."
    "O[B:18](O)[c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1"
    ">>[CH3:1][C:2]([CH3:3])([CH3:4])[c:5]1[cH:6][cH:7]"
    "[c:8](-[c:9]2[cH:10][cH:11][cH:12][cH:13][cH:14]2)"
    "[cH:15][cH:16]1"
)
KEKULIZED_SUZUKI_REACTION = (
    "COC(=O)C1=CC=C(Br)C=C1.OB(O)C1=CC=C(F)C=C1"
    ">>COC(=O)c1ccc(-c2ccc(F)cc2)cc1"
)
KEKULIZED_SUZUKI_MAPPING = (
    "[Br:18][c:8]1[cH:7][cH:6][c:5]([C:3]([O:2][CH3:1])=[O:4])"
    "[cH:17][cH:16]1.O[B:19](O)[c:9]1[cH:10][cH:11][c:12]([F:13])"
    "[cH:14][cH:15]1>>[CH3:1][O:2][C:3](=[O:4])[c:5]1[cH:6]"
    "[cH:7][c:8](-[c:9]2[cH:10][cH:11][c:12]([F:13])[cH:14]"
    "[cH:15]2)[cH:16][cH:17]1"
)
REPEATED_SUZUKI_REACTION = (
    "Brc1ccc(Br)cc1.OB(O)c1ccccc1.OB(O)c1ccccc1"
    ">>c1ccc(-c2ccc(-c3ccccc3)cc2)cc1"
)
REPEATED_SUZUKI_MAPPING = (
    "[Br:19][c:5]1[cH:6][cH:7][c:8]([Br:20])[cH:15][cH:16]1."
    "O[B:21](O)[c:4]1[cH:3][cH:2][cH:1][cH:18][cH:17]1."
    "O[B:22](O)[c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1"
    ">>[cH:1]1[cH:2][cH:3][c:4](-[c:5]2[cH:6][cH:7]"
    "[c:8](-[c:9]3[cH:10][cH:11][cH:12][cH:13][cH:14]3)"
    "[cH:15][cH:16]2)[cH:17][cH:18]1"
)
REORDERED_STEREO_REACTION = (
    "Brc1ccccn1.N[C@@H]1C[C@@H](OC)C1"
    ">>c1cccc(N[C@@H]2C[C@@H](OC)C2)n1"
)
REORDERED_STEREO_MAPPING = (
    "[Br:14][c:7]1[cH:8][cH:9][cH:10][cH:11][n:12]1."
    "[CH3:1][O:2][C@H:3]1[CH2:4][C@H:5]([NH2:6])[CH2:13]1"
    ">>[CH3:1][O:2][C@H:3]1[CH2:4][C@H:5]("
    "[NH:6][c:7]2[cH:8][cH:9][cH:10][cH:11][n:12]2)[CH2:13]1"
)
ACYL_CHLORIDE_REACTION = "CCCCCCCCC(=O)O>>CCCCCCCCC(=O)Cl"
ACYL_CHLORIDE_MAPPING = (
    "[OH:12][C:9]([CH2:8][CH2:7][CH2:6][CH2:5]"
    "[CH2:4][CH2:3][CH2:2][CH3:1])=[O:10]"
    ">>[CH3:1][CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]"
    "[CH2:7][CH2:8][C:9](=[O:10])[Cl:11]"
)
SUBSTITUTED_ACYL_CHLORIDE_REACTION = (
    "Cc1c(F)c([N+](=O)[O-])cc(C(=O)O)c1Cl"
    ">>Cc1c(F)c([N+](=O)[O-])cc(C(=O)Cl)c1Cl"
)
SUBSTITUTED_ACYL_CHLORIDE_MAPPING = (
    "[OH:16][C:11]([c:10]1[cH:9][c:5]([N+:6](=[O:7])[O-:8])"
    "[c:3]([F:4])[c:2]([CH3:1])[c:14]1[Cl:13])=[O:12]"
    ">>[CH3:1][c:2]1[c:3]([F:4])[c:5]([N+:6](=[O:7])[O-:8])"
    "[cH:9][c:10]([C:11](=[O:12])[Cl:13])[c:14]1[Cl:15]"
)


class _ResolvedFixtureProvider:
    metadata = METADATA

    def __init__(
        self,
        *,
        mapped_reaction=RESOLVED_SUZUKI_MAPPING,
        allowed_pairs=frozenset({("B", "C"), ("C", "C"), ("C", "O")}),
        signature_conflict=False,
    ):
        self.mapped_reaction = mapped_reaction
        self.allowed_pairs = allowed_pairs
        self.signature_conflict = signature_conflict

    def map_reactions(self, reactions):
        results = []
        for reaction in reactions:
            result = validate_external_atom_mapping(
                reaction,
                self.mapped_reaction,
                provider_metadata=self.metadata,
                mapper_confidence=0.44,
            )
            assert result.normalization is not None
            core_edits = tuple(
                edit
                for edit in result.normalization.edits
                if tuple(
                    sorted(
                        (
                            edit.atom_1.element,
                            edit.atom_2.element if edit.atom_2 else "H",
                        )
                    )
                )
                in self.allowed_pairs
            )
            if self.signature_conflict:
                core_edits = tuple(
                    edit for edit in core_edits if edit.edit_type != "formed"
                )
            results.append(
                replace(
                    result,
                    normalization=replace(
                        result.normalization,
                        edits=core_edits,
                    ),
                )
            )
        return tuple(results)


def test_external_mapping_preserves_structure_and_extracts_fischer_edits() -> None:
    result = validate_external_atom_mapping(
        FISCHER_REACTION,
        FISCHER_RXNMAPPER_OUTPUT,
        provider_metadata=METADATA,
        mapper_confidence=0.65,
    )

    assert result.valid
    assert result.structure_preserved
    assert result.product_mapping_coverage == 1.0
    assert result.normalization is not None
    assert result.normalization.evidence == EXTERNAL_MAPPING_EVIDENCE
    assert len(result.normalization.edits) == 9
    edit_types = [edit.edit_type for edit in result.normalization.edits]
    assert edit_types.count("formed") == 2
    assert edit_types.count("broken") == 2
    assert edit_types.count("order_changed") == 2
    assert edit_types.count("hydrogen_change") == 3
    assert all(
        edit.evidence == EXTERNAL_MAPPING_EVIDENCE
        and edit.confidence == 0.65
        for edit in result.normalization.edits
    )
    assert (
        "EXTERNAL_MAPPING_PROJECTED_REACTANT_BOUNDARY_ATOMS:2"
        in result.warnings
    )


def test_boundary_projection_does_not_map_unattached_spectator() -> None:
    result = validate_external_atom_mapping(
        FISCHER_REACTION,
        FISCHER_RXNMAPPER_OUTPUT,
        provider_metadata=METADATA,
        mapper_confidence=0.65,
    )

    assert result.mapped_reaction_smiles is not None
    parsed = parse_reaction_smiles(result.mapped_reaction_smiles)
    chloride = next(
        component
        for component in parsed.reactants
        if component.canonical_smiles == "Cl"
    )
    molecule = parse_smiles(chloride.input_smiles)
    assert molecule is not None
    assert molecule.GetAtomWithIdx(0).GetAtomMapNum() == 0


def test_external_mapping_rejects_changed_chemistry() -> None:
    result = validate_external_atom_mapping(
        "CCO>>CC=O",
        "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][OH:3]",
        provider_metadata=METADATA,
        mapper_confidence=0.99,
    )

    assert not result.valid
    assert not result.structure_preserved
    assert result.normalization is None
    assert result.error == "MAPPER_CHANGED_REACTION_STRUCTURE"
    assert "EXTERNAL_MAPPER_CHANGED_REACTION_STRUCTURE" in result.warnings


def test_external_mapping_normalizes_stereo_after_atom_reordering() -> None:
    result = validate_external_atom_mapping(
        REORDERED_STEREO_REACTION,
        REORDERED_STEREO_MAPPING,
        provider_metadata=METADATA,
        mapper_confidence=0.8,
    )

    assert result.valid
    assert result.structure_preserved
    assert result.normalization is not None
    assert {
        (edit.edit_type, edit.atom_1.element, edit.atom_2.element)
        for edit in result.normalization.edits
        if edit.atom_2 is not None
    } == {
        ("broken", "C", "Br"),
        ("formed", "N", "C"),
    }


def test_external_mapping_still_rejects_a_real_stereo_change() -> None:
    result = validate_external_atom_mapping(
        "F[C@H](Cl)Br>>F[C@H](Cl)I",
        "[F:1][C@@H:2]([Cl:3])[Br:4]>>[F:1][C@@H:2]([Cl:3])[I:4]",
        provider_metadata=METADATA,
        mapper_confidence=0.8,
    )

    assert not result.valid
    assert not result.structure_preserved
    assert result.error == "MAPPER_CHANGED_REACTION_STRUCTURE"


def test_signature_unavailable_core_uses_original_acyl_coordinates() -> None:
    base = featurize_reaction(ACYL_CHLORIDE_REACTION)

    assessment = analyze_reaction_with_external_mapping(
        ACYL_CHLORIDE_REACTION,
        _ResolvedFixtureProvider(
            mapped_reaction=ACYL_CHLORIDE_MAPPING,
            allowed_pairs=frozenset({("C", "Cl"), ("C", "O")}),
        ),
        base_analysis=base,
    )

    assert assessment.status == "external_mapping_signature_unavailable"
    assert assessment.analysis.reaction_core is not None
    assert {
        (
            transition.atom_map_number,
            (
                transition.before_state.atom_index
                if transition.before_state is not None
                else None
            ),
            (
                transition.after_state.atom_index
                if transition.after_state is not None
                else None
            ),
        )
        for transition in assessment.analysis.reaction_core.atom_transitions
    } == {
        (9, 8, 8),
        (11, None, 10),
        (12, 10, None),
    }

    projection = build_reaction_display_projection(
        reaction_render_context_from_analysis(assessment.analysis)
    )
    assert projection.minimum_reaction_smiles == (
        "*C(=O)O>>*C(=O)Cl"
    )
    assert projection.render_reaction_smiles == (
        "O=C(O)[*:1]>>O=C(Cl)[*:1]"
    )


def test_external_mapping_repairs_unique_spectator_chlorine_swap() -> None:
    mapping = validate_external_atom_mapping(
        SUBSTITUTED_ACYL_CHLORIDE_REACTION,
        SUBSTITUTED_ACYL_CHLORIDE_MAPPING,
        provider_metadata=METADATA,
        mapper_confidence=0.68,
    )

    assert mapping.valid
    assert mapping.normalization is not None
    assert "EXTERNAL_MAPPING_RECONCILED_EQUIVALENT_ATOMS:2" in (
        mapping.warnings
    )
    assert {
        (
            edit.edit_type,
            tuple(
                sorted(
                    (
                        edit.atom_1.element,
                        edit.atom_2.element if edit.atom_2 else "H",
                    )
                )
            ),
        )
        for edit in mapping.normalization.edits
    } == {
        ("broken", ("C", "O")),
        ("formed", ("C", "Cl")),
    }

    assessment = analyze_reaction_with_external_mapping(
        SUBSTITUTED_ACYL_CHLORIDE_REACTION,
        _ResolvedFixtureProvider(
            mapped_reaction=SUBSTITUTED_ACYL_CHLORIDE_MAPPING,
            allowed_pairs=frozenset({("C", "Cl"), ("C", "O")}),
        ),
        base_analysis=featurize_reaction(
            SUBSTITUTED_ACYL_CHLORIDE_REACTION
        ),
    )
    assert assessment.analysis.reaction_core is not None
    assert "REACTION_CORE_NO_OP_PRIMARY_CENTER" not in (
        assessment.analysis.reaction_core.warnings
    )
    projection = build_reaction_display_projection(
        reaction_render_context_from_analysis(assessment.analysis)
    )
    assert projection.minimum_reaction_smiles == (
        "*C(=O)O>>*C(=O)Cl"
    )
    assert projection.render_reaction_smiles == (
        "O=C(O)[*:1]>>O=C(Cl)[*:1]"
    )


def test_external_mapping_preserves_ambiguous_chlorine_reassignment() -> None:
    reaction = "CC(Cl)Cl>>ClCC(Cl)Cl"
    ambiguous_mapping = (
        "[CH3:4][CH:6]([Cl:1])[Cl:2]"
        ">>[CH2:4]([Cl:2])[CH:6]([Cl:3])[Cl:1]"
    )

    result = validate_external_atom_mapping(
        reaction,
        ambiguous_mapping,
        provider_metadata=METADATA,
        mapper_confidence=0.5,
    )

    assert result.valid
    assert result.mapped_reaction_smiles == ambiguous_mapping
    assert "EXTERNAL_MAPPING_RECONCILIATION_AMBIGUOUS:Cl:2" in (
        result.warnings
    )
    assert not any(
        warning.startswith("EXTERNAL_MAPPING_RECONCILED_EQUIVALENT_ATOMS")
        for warning in result.warnings
    )


def test_external_mapping_does_not_rewrite_a_minimal_assignment() -> None:
    reaction = "CC(Cl)Cl>>ClCC(Cl)Cl"
    minimal_mapping = (
        "[CH3:4][CH:6]([Cl:1])[Cl:2]"
        ">>[CH2:4]([Cl:3])[CH:6]([Cl:1])[Cl:2]"
    )

    result = validate_external_atom_mapping(
        reaction,
        minimal_mapping,
        provider_metadata=METADATA,
        mapper_confidence=0.5,
    )

    assert result.valid
    assert result.mapped_reaction_smiles == minimal_mapping
    assert not any(
        "RECONCIL" in warning for warning in result.warnings
    )


class _FakeMapper:
    def get_attention_guided_atom_maps(self, reactions: list[str]) -> list[dict]:
        return [
            {
                "mapped_rxn": FISCHER_RXNMAPPER_OUTPUT,
                "confidence": 0.65,
            }
            for _ in reactions
        ]


def test_rxnmapper_provider_construction_does_not_import_mapper(monkeypatch) -> None:
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name == "rxnmapper" or name.startswith("rxnmapper."):
            raise AssertionError("provider construction imported RXNMapper")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    provider = RxnMapperProvider(mapper=_FakeMapper())

    assert provider.metadata.provider_id == "rxnmapper"


def test_rxnmapper_provider_batches_and_retains_model_provenance() -> None:
    provider = RxnMapperProvider(batch_size=1, mapper=_FakeMapper())
    results = provider.map_reactions([FISCHER_REACTION, FISCHER_REACTION])

    assert len(results) == 2
    assert all(result.valid for result in results)
    assert provider.metadata.provider_id == "rxnmapper"
    assert provider.metadata.provider_version
    assert provider.metadata.model_id == "albert_heads_8_uspto_all_1310k"


def test_forced_resolved_mapping_only_enriches_the_minimized_core() -> None:
    base = featurize_reaction(RESOLVED_SUZUKI_REACTION)

    assessment = analyze_reaction_with_external_mapping(
        RESOLVED_SUZUKI_REACTION,
        _ResolvedFixtureProvider(),
        base_analysis=base,
        force_resolved_shadow=True,
    )

    assert assessment.status == "external_mapping_internal_consensus"
    assert assessment.analysis.reaction_label.text
    assert "EXTERNAL_MAPPING_INTERNAL_CONSENSUS" in (
        assessment.analysis.reaction_label.warnings
    )
    assert assessment.analysis.reaction_signature == base.reaction_signature
    assert assessment.analysis.reaction_core is not None
    assert assessment.analysis.reaction_core.evidence_status == "external"
    assert assessment.analysis.reaction_core.evidence == (
        "external_mapping_internal_consensus"
    )
    assert not hasattr(assessment.analysis.reaction_core, "generic_label")
    for transition in assessment.analysis.reaction_core.atom_transitions:
        for state, components in (
            (transition.before_state, base.reactants),
            (transition.after_state, base.products),
        ):
            if state is None:
                continue
            molecule = parse_smiles(
                components[state.component_index].input_smiles
            )
            assert molecule is not None
            assert molecule.GetAtomWithIdx(state.atom_index).GetSymbol() == (
                state.element
            )
    assert "REACTION_CORE_EXTERNAL_MAPPING_PROPOSAL" in (
        assessment.analysis.reaction_core.warnings
    )
    assert build_reaction_review_summary(
        assessment.analysis
    ).reaction_label


def test_forced_resolved_mapping_rerenders_material_core_context() -> None:
    reaction = "CSC.[O]>>CS(C)=O"
    mapped = (
        "[CH3:1][S:2][CH3:3].[O:4]>>"
        "[CH3:1][S:2]([CH3:3])=[O:4]"
    )
    base = featurize_reaction(reaction)

    assessment = analyze_reaction_with_external_mapping(
        reaction,
        _ResolvedFixtureProvider(
            mapped_reaction=mapped,
            allowed_pairs=frozenset({("O", "S")}),
        ),
        base_analysis=base,
        force_resolved_shadow=True,
    )

    assert assessment.status == "external_mapping_internal_consensus"
    assert assessment.analysis.reaction_signature == base.reaction_signature
    assert assessment.analysis.reaction_core is not None
    assert assessment.analysis.observation is not None
    assert (
        assessment.analysis.observation.core
        == assessment.analysis.reaction_core
    )
    assert assessment.analysis.reaction_label.text == "O + S → S–Alk–Alk=O"
    assert assessment.analysis.reaction_label.basis == "reaction_sites"
    assert assessment.analysis.reaction_label.evidence == base.evidence_quality


def test_forced_mapping_projects_maps_without_smiles_coordinate_round_trip() -> None:
    base = featurize_reaction(KEKULIZED_SUZUKI_REACTION)
    assert base.reaction_signature is not None
    assert base.reaction_core is not None

    assessment = analyze_reaction_with_external_mapping(
        KEKULIZED_SUZUKI_REACTION,
        _ResolvedFixtureProvider(
            mapped_reaction=KEKULIZED_SUZUKI_MAPPING,
            allowed_pairs=frozenset(
                {("B", "C"), ("Br", "C"), ("C", "C")}
            ),
        ),
        base_analysis=base,
        force_resolved_shadow=True,
    )

    assert assessment.status in {
        "external_mapping_internal_consensus",
        "external_mapping_only",
    }
    assert assessment.analysis.reaction_signature == base.reaction_signature
    assert assessment.analysis.reaction_core is not None
    assert assessment.analysis.reaction_core.evidence_status == "external"
    assert "EXTERNAL_MAPPING_CORE_COORDINATE_PROJECTION_FAILED" not in (
        assessment.analysis.warnings
    )


def test_forced_mapping_conflict_retains_the_resolved_interpretation() -> None:
    base = featurize_reaction(RESOLVED_SUZUKI_REACTION)

    assessment = analyze_reaction_with_external_mapping(
        RESOLVED_SUZUKI_REACTION,
        _ResolvedFixtureProvider(signature_conflict=True),
        base_analysis=base,
        force_resolved_shadow=True,
    )

    assert assessment.status == "external_mapping_signature_conflict"
    assert assessment.analysis.reaction_label == base.reaction_label
    assert assessment.analysis.reaction_label == base.reaction_label
    assert assessment.analysis.reaction_signature == base.reaction_signature
    assert assessment.analysis.reaction_core is not None
    assert assessment.analysis.reaction_core.evidence_status == "hypothesis"
    assert "REACTION_CORE_CONFLICTING_EVIDENCE_PROPOSAL" in (
        assessment.analysis.reaction_core.warnings
    )
    assert "EXTERNAL_MAPPING_SIGNATURE_CONFLICT" in assessment.analysis.warnings


def test_forced_mapping_builds_one_core_with_two_suzuki_events() -> None:
    base = featurize_reaction(REPEATED_SUZUKI_REACTION)

    assessment = analyze_reaction_with_external_mapping(
        REPEATED_SUZUKI_REACTION,
        _ResolvedFixtureProvider(
            mapped_reaction=REPEATED_SUZUKI_MAPPING,
            allowed_pairs=frozenset(
                {("B", "C"), ("Br", "C"), ("C", "C")}
            ),
        ),
        base_analysis=base,
        force_resolved_shadow=True,
    )

    assert assessment.status in {
        "external_mapping_internal_consensus",
        "external_mapping_only",
    }
    assert base.reaction_signature is None
    assert assessment.analysis.reaction_signature is not None
    assert assessment.analysis.reaction_label.text
    assert assessment.analysis.reaction_label.status == "multi_event"
    assert assessment.analysis.reaction_core is not None
    assert assessment.analysis.reaction_core.event_count == 2
    assert not hasattr(assessment.analysis.reaction_core, "generic_label")
