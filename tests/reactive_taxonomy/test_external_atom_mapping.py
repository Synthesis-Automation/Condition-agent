"""Tests for optional external atom-mapping validation and provenance."""

from __future__ import annotations

from dataclasses import replace

from reactive_taxonomy import (
    AtomMappingProviderMetadata,
    EXTERNAL_MAPPING_EVIDENCE,
    RxnMapperProvider,
    analyze_reaction_with_external_mapping,
    build_reaction_review_summary,
    featurize_reaction,
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


class _FakeMapper:
    def get_attention_guided_atom_maps(self, reactions: list[str]) -> list[dict]:
        return [
            {
                "mapped_rxn": FISCHER_RXNMAPPER_OUTPUT,
                "confidence": 0.65,
            }
            for _ in reactions
        ]


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
    assert assessment.analysis.reaction_label == base.reaction_label
    assert assessment.analysis.reaction_label_status == base.reaction_label_status
    assert assessment.analysis.display_label == base.display_label
    assert assessment.analysis.reaction_signature == base.reaction_signature
    assert assessment.analysis.reaction_core is not None
    assert assessment.analysis.reaction_core.evidence_status == "external"
    assert assessment.analysis.reaction_core.evidence == (
        "external_mapping_internal_consensus"
    )
    assert assessment.analysis.reaction_core.generic_label == (
        "Ar–B + Ar–O → Ar–Ar"
    )
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
    ).detailed_reaction_label == build_reaction_review_summary(
        base
    ).detailed_reaction_label


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
    assert assessment.analysis.display_label == base.display_label
    assert assessment.analysis.reaction_signature == base.reaction_signature
    assert assessment.analysis.reaction_core is None
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

    assert assessment.status == "external_mapping_internal_consensus"
    assert assessment.analysis.reaction_signature == base.reaction_signature
    assert assessment.analysis.display_label == base.display_label
    assert assessment.analysis.reaction_core is not None
    assert assessment.analysis.reaction_core.event_count == 2
    assert assessment.analysis.reaction_core.generic_label == (
        "2 × Ar–B + 2 × Ar–Br → 2 × Ar–Ar"
    )
