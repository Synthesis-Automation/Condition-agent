"""Converter and recommender integration tests for optional atom mapping."""

from __future__ import annotations

import json
from dataclasses import dataclass

from reactive_taxonomy import (
    AtomMappingProviderMetadata,
    ExternalAtomMappingResult,
    analyze_reaction_with_external_mapping,
    validate_external_atom_mapping,
)

from condition_recommender import GenericConditionRecommender
from condition_recommender.conversion.engine import convert_datasets
from condition_recommender.conversion.concise_review import (
    concise_reaction_review_row,
)
from condition_recommender.conversion.generic import (
    GenericConversionCache,
    convert_record,
)
import condition_recommender.conversion.sharded as sharded_module
from condition_recommender.conversion.sharded import convert_datasets_sharded
from condition_recommender.conversion.input_schema import adapt_row
from condition_recommender.generic_indexing import build_generic_index
from condition_recommender.models import (
    AdmissionTier,
    ChemistryStatus,
    IndexEligibility,
)


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
CONFLICT_REACTION = "Cl.NNc1ccccc1.O=C1CCCCC1>>c1ccc2c3c([nH]c2c1)CCCC3"
CONFLICT_RXNMAPPER_OUTPUT = (
    "Cl.[NH2:14][NH:7][c:8]1[cH:4][cH:3][cH:2][cH:1][cH:9]1."
    "[O:15]=[C:6]1[CH2:5][CH2:10][CH2:11][CH2:12][CH2:13]1"
    ">>[cH:1]1[cH:2][cH:3][c:4]2[c:5]3[c:6]"
    "([nH:7][c:8]2[cH:9]1)[CH2:10][CH2:11][CH2:12][CH2:13]3"
)
MAPPER_ONLY_REACTION = (
    "CI.CC(=NNc1ccccc1)c1ccccc1"
    ">>Cn1c(-c2ccccc2)cc2ccccc21"
)
MAPPER_ONLY_RXNMAPPER_OUTPUT = (
    "[NH:17]([N:2]=[C:3]([c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1)"
    "[CH3:10])[c:11]1[cH:12][cH:13][cH:14][cH:15][cH:16]1."
    "[I:18][CH3:1]>>[CH3:1][n:2]1[c:3](-[c:4]2[cH:5][cH:6]"
    "[cH:7][cH:8][cH:9]2)[cH:10][c:11]2[cH:12][cH:13][cH:14]"
    "[cH:15][c:16]12"
)
INCOMPLETE_ACETAL_REACTION = (
    "CC(C)=O.OCC(O)C(O)C(O)CO"
    ">>CC1(C)OC[C@@H]([C@H]2OC(C)(C)O[C@H]2CO)O1"
)
INCOMPLETE_ACETAL_RXNMAPPER_OUTPUT = (
    "[CH3:1][C:9](=[O:8])[CH3:10]."
    "[OH:17][CH:7]([CH:6]([CH2:5][OH:4])[OH:16])"
    "[CH:13]([OH:12])[CH2:14][OH:15]"
    ">>[CH3:1][C:2]1([CH3:3])[O:4][CH2:5][C@@H:6]"
    "([C@H:7]2[O:8][C:9]([CH3:10])([CH3:11])[O:12]"
    "[C@H:13]2[CH2:14][OH:15])[O:16]1"
)
ACYL_FLUORIDE_REACTION = "O=C(O)c1ccccc1>>O=C(F)c1ccccc1"
ACYL_FLUORIDE_RXNMAPPER_OUTPUT = (
    "[OH:10][C:2](=[O:1])[c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1"
    ">>[O:1]=[C:2]([F:3])[c:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1"
)
METADATA = AtomMappingProviderMetadata(
    provider_id="fixture_mapper",
    provider_version="1.0",
    model_id="fixture",
    model_sha256="abc",
)


@dataclass
class _FixtureProvider:
    mapped_reaction: str = FISCHER_RXNMAPPER_OUTPUT
    confidence: float = 0.65
    call_count: int = 0
    metadata: AtomMappingProviderMetadata = METADATA

    def map_reactions(
        self,
        reaction_smiles: tuple[str, ...],
    ) -> tuple[ExternalAtomMappingResult, ...]:
        self.call_count += 1
        return tuple(
            validate_external_atom_mapping(
                reaction,
                self.mapped_reaction,
                provider_metadata=self.metadata,
                mapper_confidence=self.confidence,
            )
            for reaction in reaction_smiles
        )


class _ShardedFixtureProvider(_FixtureProvider):
    @staticmethod
    def is_available() -> bool:
        return True


def _raw_record():
    return adapt_row(
        {
            "reaction_id": "fischer-1",
            "reaction_type": "Fischer indole synthesis",
            "reaction_smiles": FISCHER_REACTION,
            "yield_pct": "82",
            "reagent_cas": "7647-01-0",
            "solvent_cas": "64-17-5",
            "reference": "Fixture reference",
            "stages": "1",
        },
        source_dataset="fischer",
        source_path="fischer.csv",
        source_row_number=2,
    )


def test_external_mapping_selects_one_internal_fischer_hypothesis() -> None:
    assessment = analyze_reaction_with_external_mapping(
        FISCHER_REACTION,
        _FixtureProvider(),
    )

    assert assessment.status == "external_mapping_internal_consensus"
    assert len(assessment.matched_hypothesis_ids) == 1
    assert assessment.analysis.reaction_signature is not None
    assert (
        assessment.analysis.evidence_quality
        == "external_mapping_internal_consensus"
    )
    assert assessment.analysis.named_family is None
    assert assessment.analysis.edit_hypotheses


def test_external_mapping_conflict_retains_internal_hypotheses() -> None:
    assessment = analyze_reaction_with_external_mapping(
        CONFLICT_REACTION,
        _FixtureProvider(
            mapped_reaction=CONFLICT_RXNMAPPER_OUTPUT,
            confidence=0.31,
        ),
    )

    assert assessment.status == "external_mapping_hypothesis_conflict"
    assert assessment.analysis.reaction_signature is None
    assert len(assessment.analysis.edit_hypotheses) == 2
    assert "EXTERNAL_MAPPING_HYPOTHESIS_CONFLICT" in assessment.analysis.warnings


def test_mapper_only_edits_create_review_qualified_signature() -> None:
    assessment = analyze_reaction_with_external_mapping(
        MAPPER_ONLY_REACTION,
        _FixtureProvider(
            mapped_reaction=MAPPER_ONLY_RXNMAPPER_OUTPUT,
            confidence=0.46,
        ),
    )

    assert assessment.status == "external_mapping_only"
    assert assessment.matched_hypothesis_ids == ()
    assert assessment.analysis.reaction_signature is not None
    assert assessment.analysis.evidence_quality == "external_atom_mapping"
    assert "EXTERNAL_MAPPING_WITHOUT_INTERNAL_CONSENSUS" in (
        assessment.analysis.warnings
    )


def test_product_only_mapped_endpoints_do_not_break_reaction_topology() -> None:
    assessment = analyze_reaction_with_external_mapping(
        INCOMPLETE_ACETAL_REACTION,
        _FixtureProvider(
            mapped_reaction=INCOMPLETE_ACETAL_RXNMAPPER_OUTPUT,
            confidence=0.35,
        ),
    )

    assert assessment.status == "external_mapping_signature_unavailable"
    assert assessment.analysis.valid
    assert assessment.analysis.reaction_signature is None
    assert assessment.analysis.reaction_completeness is not None
    assert assessment.analysis.reaction_completeness.status == "incomplete"
    assert "EXTERNAL_MAPPING_SIGNATURE_UNAVAILABLE" in assessment.warnings


def test_mapper_unavailable_signature_does_not_block_supported_partial_record() -> None:
    raw = adapt_row(
        {
            "reaction_id": "acyl-fluoride-1",
            "reaction_smiles": ACYL_FLUORIDE_REACTION,
            "yield_pct": "90",
            "reagent_cas": "63517-29-3",
            "solvent_cas": "141-78-6",
            "reference": "Acyl fluoride reference",
        },
        source_dataset="acyl-fluoride",
        source_path="acyl-fluoride.csv",
        source_row_number=1,
    )
    record = convert_record(
        raw,
        mapping_provider=_FixtureProvider(
            mapped_reaction=ACYL_FLUORIDE_RXNMAPPER_OUTPUT,
            confidence=0.88,
        ),
    )

    assert record.external_atom_mapping is not None
    assert record.external_atom_mapping["status"] == (
        "external_mapping_signature_unavailable"
    )
    assert record.chemistry_status == ChemistryStatus.REVIEW
    assert record.index_eligibility == IndexEligibility.ELIGIBLE
    assert record.fragment_source_support[0]["status"] == "supported"
    assert len(build_generic_index([record.to_dict()]).rows) == 1


def test_converter_persists_mapper_provenance_but_excludes_precedent() -> None:
    provider = _FixtureProvider()
    record = convert_record(_raw_record(), mapping_provider=provider)

    assert record.reaction_signature is not None
    assert record.external_atom_mapping is not None
    assert (
        record.external_atom_mapping["status"]
        == "external_mapping_internal_consensus"
    )
    assert record.external_atom_mapping["provider"]["provider_id"] == (
        "fixture_mapper"
    )
    assert record.chemistry_status == ChemistryStatus.REVIEW
    assert record.admission_tier == AdmissionTier.REVIEW
    assert record.index_eligibility == IndexEligibility.REVIEW_ONLY
    assert "external_mapping_review_required" in record.admission_reasons
    review = concise_reaction_review_row(record.to_dict())
    assert review["external_mapping_status"] == (
        "external_mapping_internal_consensus"
    )
    assert review["external_mapping_provider"] == "fixture_mapper"
    assert review["external_mapping_confidence"] == "0.65"
    assert not build_generic_index([record.to_dict()]).rows
    incorrectly_promoted = record.to_dict()
    incorrectly_promoted["index_eligibility"] = "eligible"
    assert not build_generic_index([incorrectly_promoted]).rows


def test_converter_cache_maps_duplicate_reactions_once() -> None:
    provider = _FixtureProvider()
    cache = GenericConversionCache()

    first = convert_record(
        _raw_record(),
        cache=cache,
        mapping_provider=provider,
    )
    second = convert_record(
        _raw_record(),
        cache=cache,
        mapping_provider=provider,
    )

    assert provider.call_count == 1
    assert first.reaction_signature == second.reaction_signature
    assert first.external_atom_mapping == second.external_atom_mapping


def test_conversion_engine_reports_external_mapping_dispositions(tmp_path) -> None:
    source = tmp_path / "fischer.csv"
    source.write_text(
        "reaction_id,reaction_type,reaction_smiles,yield_pct,reagent_cas,"
        "solvent_cas,reference,stages\n"
        f'fischer-1,Fischer indole synthesis,"{FISCHER_REACTION}",82,'
        "7647-01-0,64-17-5,Fixture reference,1\n",
        encoding="utf-8",
    )
    provider = _FixtureProvider()

    report = convert_datasets(
        source,
        tmp_path / "converted",
        mapping_provider=provider,
    )

    assert report["schema_version"] == "2.1"
    assert report["signature_count"] == 1
    assert report["index_eligibility_counts"]["review_only"] == 1
    assert report["external_atom_mapping"] == {
        "enabled": True,
        "provider": {
            "provider_id": "fixture_mapper",
            "provider_version": "1.0",
            "model_id": "fixture",
            "model_sha256": "abc",
        },
        "status_counts": {"external_mapping_internal_consensus": 1},
    }
    review_header = (tmp_path / "converted" / "review.csv").read_text(
        encoding="utf-8-sig",
    ).splitlines()[0]
    assert "external_mapping_status" in review_header
    assert "external_atom_mapping_json" in review_header


def test_sharded_converter_carries_mapping_into_manifest_and_records(
    tmp_path,
    monkeypatch,
) -> None:
    source = tmp_path / "fischer.csv"
    source.write_text(
        "reaction_id,reaction_type,reaction_smiles,yield_pct,reagent_cas,"
        "solvent_cas,reference,stages\n"
        f'fischer-1,Fischer indole synthesis,"{FISCHER_REACTION}",82,'
        "7647-01-0,64-17-5,Fixture reference,1\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        sharded_module,
        "RxnMapperProvider",
        _ShardedFixtureProvider,
    )

    report = convert_datasets_sharded(
        source,
        tmp_path / "sharded",
        shard_size=1,
        workers=1,
        use_rxnmapper=True,
    )

    assert report["external_atom_mapping"]["enabled"]
    assert report["external_atom_mapping"]["provider_id"] == "fixture_mapper"
    assert report["external_atom_mapping"]["status_counts"] == {
        "external_mapping_internal_consensus": 1
    }
    assert report["index_eligibility_counts"] == {"review_only": 1}
    manifest = json.loads(
        (tmp_path / "sharded" / "shard_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["definition_contract"]["external_atom_mapping"] == {
        "enabled": True,
        "provider_id": "fixture_mapper",
        "provider_version": "1.0",
        "model_id": "fixture",
        "model_sha256": "abc",
    }
    assert report["integrity"]["valid"]


def test_recommender_uses_mapper_supported_query_with_review_cautions() -> None:
    converted = convert_record(
        _raw_record(),
        mapping_provider=_FixtureProvider(),
    ).to_dict()
    converted.update(
        {
            "reaction_id": "verified-precedent",
            "observation_id": "verified-observation",
            "index_eligibility": "eligible",
            "chemistry_status": "verified",
            "admission_tier": "verified",
            "external_atom_mapping": None,
        }
    )
    index = build_generic_index([converted])
    provider = _FixtureProvider()

    result = GenericConditionRecommender(
        index=index,
        mapping_provider=provider,
    ).recommend(FISCHER_REACTION, minimum_pool_size=1)

    assert result.valid
    assert result.external_mapping_status == (
        "external_mapping_internal_consensus"
    )
    assert result.external_mapping_provider == "fixture_mapper"
    assert result.external_mapping_confidence == 0.65
    assert result.recommendation_mode == "external_mapping_internal_consensus"
    assert "QUERY_SIGNATURE_USES_EXTERNAL_MAPPING" in result.warnings
    assert result.recommendations
    assert any(
        "external mapper" in caution
        for caution in result.recommendations[0].cautions
    )
