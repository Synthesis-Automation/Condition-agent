"""Preserve directed molecular bonds across conversion and source-linked search."""

import json

from condition_recommender.conversion.identities import canonical_reaction_identity
from condition_recommender.fragment_index import build_fragment_index
from condition_recommender.fragment_search import lookup_exact_product
from reactive_taxonomy.chemistry.reaction_smiles import split_reaction_fields
from reactive_taxonomy.reaction_parser import parse_reaction_smiles


def test_coordinate_bonds_survive_identity_and_product_retrieval(tmp_path):
    reaction = "N->[Cu]>O>N->[Cu]Cl"
    assert split_reaction_fields(reaction) == ["N->[Cu]", "O", "N->[Cu]Cl"]
    assert split_reaction_fields("N->[Cu]>>CO") == ["N->[Cu]", "", "CO"]
    parsed = parse_reaction_smiles(reaction, include_molecular_interpretation=False)
    assert parsed.valid and len(parsed.agents) == 1
    canonical = canonical_reaction_identity(reaction)
    reversed_bond_text = canonical_reaction_identity("[Cu]<-N>O>Cl[Cu]<-N")
    assert canonical is not None and canonical == reversed_bond_text

    source = tmp_path / "records.jsonl"
    source.write_text(json.dumps({"observation_id": "coordinate", "reaction_id": canonical.reaction_id,
                                  "reaction_smiles": reaction}) + "\n", encoding="utf-8")
    index = tmp_path / "index.sqlite"
    built = build_fragment_index(source, index)
    assert built["counts"]["indexed_product_components"] == 1
    match = lookup_exact_product(index, "Cl[Cu]<-N")
    assert match["precedents"][0]["observation_id"] == "coordinate"


def test_invalid_field_count_still_rejected_with_coordinate_bond():
    reaction = "N->[Cu]>O>CO>O"
    assert not parse_reaction_smiles(reaction).valid
    assert canonical_reaction_identity(reaction) is None
