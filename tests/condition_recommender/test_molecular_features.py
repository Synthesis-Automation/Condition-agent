"""Observation-aligned molecular features through conversion, storage and ranking."""

from copy import deepcopy
from dataclasses import asdict, replace
from pathlib import Path

import pytest

from reactive_taxonomy import featurize_reaction
from condition_recommender.conversion.generic import convert_record
from condition_recommender.conversion.input_schema import adapt_row
from condition_recommender.generic_api import recommend_indexed_signature
from condition_recommender.generic_indexing import build_generic_index, build_generic_index_from_rows
from condition_recommender.molecular_features import build_reaction_molecular_features, validate_molecular_features
from condition_recommender.signature_features import environment_profile_similarity, environment_tokens
from condition_recommender.sqlite_indexing import save_sqlite_generic_index, load_sqlite_generic_index


REACTION = "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]"


def _record():
    raw = adapt_row({
        "reaction_id": "feature-review", "reaction_smiles": REACTION,
        "reagent_cas": "584-08-7", "solvent_cas": "64-17-5",
        "yield_pct": "80", "reference": "feature-review-reference",
    }, source_dataset="development", source_path="development.csv", source_row_number=2)
    return convert_record(raw)


def test_projection_has_edit_provenance_and_does_not_mutate_signature() -> None:
    analysis = featurize_reaction(REACTION)
    before = asdict(analysis.reaction_signature)
    features = build_reaction_molecular_features(analysis).to_dict()
    validate_molecular_features(features, before)
    assert len(features["partners"]) == 2
    assert all(partner["edit_indices"] and partner["active_atom_indices"] for partner in features["partners"])
    assert environment_tokens(features)
    assert not environment_tokens(before)
    assert asdict(analysis.reaction_signature) == before


def test_ambiguous_edit_hypotheses_cannot_create_observed_features() -> None:
    analysis = featurize_reaction("O=C1CCCCC1.Cl.NNc1ccc(F)cc1>>Fc1ccc2[nH]c3c(c2c1)CCCC3")
    assert analysis.edit_hypotheses
    features = build_reaction_molecular_features(analysis)
    assert features.status == "unavailable" and not features.partners


def test_projection_rejects_conflicting_locus_and_stale_definitions() -> None:
    analysis = featurize_reaction(REACTION)
    features = build_reaction_molecular_features(analysis).to_dict()
    conflict = deepcopy(features)
    conflict["partners"][0]["active_atom_indices"] = [999]
    with pytest.raises(ValueError, match="conflicts"):
        validate_molecular_features(conflict, asdict(analysis.reaction_signature))
    stale = deepcopy(features)
    stale["definition_versions"] = {}
    with pytest.raises(ValueError, match="definitions"):
        validate_molecular_features(stale)


def test_profiles_survive_conversion_sqlite_and_live_scoring(tmp_path: Path) -> None:
    record = _record()
    index = build_generic_index([record.to_dict()])
    assert index.environment_features
    path = tmp_path / "index.sqlite"
    save_sqlite_generic_index(index, path)
    restored = load_sqlite_generic_index(path)
    row = restored.rows[0]
    assert environment_tokens(row.molecular_features) == environment_tokens(record.molecular_features)
    result = recommend_indexed_signature(
        row.signature, restored, molecular_features=row.molecular_features,
        reaction_core=row.reaction_core, query_reaction_smiles=REACTION,
        minimum_pool_size=1, top_k=1,
    )
    assert result.valid and result.recommendations
    assert result.reaction_partners[0]["reactivity_profile"]
    assert result.recommendations[0].score_trace.similarity_components["environment"] == pytest.approx(1)


def test_unresolved_profiles_never_receive_similarity_credit() -> None:
    analysis = featurize_reaction("[CH3:1][C:2]#[N:3]>>[CH3:1][CH2:2][NH2:3]")
    features = build_reaction_molecular_features(analysis).to_dict()
    assert features["partners"]
    assert all(p["reactivity_profile"]["status"] == "unresolved" for p in features["partners"])
    assert environment_profile_similarity(features, features) == 0


def test_feature_alignment_is_symmetric_and_partner_order_invariant() -> None:
    features = build_reaction_molecular_features(featurize_reaction(REACTION)).to_dict()
    reversed_features = {**features, "partners": tuple(reversed(features["partners"]))}
    assert environment_profile_similarity(features, reversed_features) == pytest.approx(1)
    assert environment_profile_similarity(reversed_features, features) == pytest.approx(1)
    row = build_generic_index([_record().to_dict()]).rows[0]
    stale = {**row.molecular_features, "definition_versions": {}}
    with pytest.raises(ValueError, match="definitions"):
        build_generic_index_from_rows([replace(row, molecular_features=stale)])
