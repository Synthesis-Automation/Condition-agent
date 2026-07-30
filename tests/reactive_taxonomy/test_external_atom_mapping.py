"""Tests for optional external atom-mapping validation and provenance."""

from __future__ import annotations

from reactive_taxonomy import (
    AtomMappingProviderMetadata,
    EXTERNAL_MAPPING_EVIDENCE,
    RxnMapperProvider,
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
