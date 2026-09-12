"""Observation refresh reuses the canonical converter and preserves source identity."""

from copy import deepcopy

import pytest

from condition_recommender.conversion.generic import convert_record
from condition_recommender.conversion.input_schema import adapt_row
from condition_recommender.conversion.refresh import (
    refresh_canonical_record,
    raw_record_from_canonical,
)


def _record():
    return convert_record(
        adapt_row(
            {
                "reaction_id": "refresh-example",
                "reaction_smiles": "[CH3:1][Br:2].[NH2:3][CH3:4]>>[CH3:1][NH:3][CH3:4]",
                "reagent_cas": "584-08-7",
                "solvent_cas": "64-17-5",
                "yield_pct": "80",
                "reference": "development-reference",
            },
            source_dataset="development",
            source_path="development.csv",
            source_row_number=2,
        )
    ).to_dict()


def test_refresh_matches_fresh_conversion_and_preserves_frozen_input():
    source = _record()
    original = deepcopy(source)
    refreshed = refresh_canonical_record(source)
    assert source == original
    assert refreshed == convert_record(raw_record_from_canonical(source)).to_dict()
    assert refreshed["reaction_signature"] == source["reaction_signature"]
    assert refreshed["observation_id"] == source["observation_id"]


def test_refresh_rejects_stale_structural_contract():
    source = _record()
    source["reaction_signature"]["definition_versions"] = {}
    with pytest.raises(ValueError, match="Structural definitions"):
        refresh_canonical_record(source)


def test_refresh_does_not_silently_erase_external_mapping_provenance():
    source = _record()
    source["external_atom_mapping"] = {"provider": "external"}
    with pytest.raises(ValueError, match="mapping provider"):
        refresh_canonical_record(source)
