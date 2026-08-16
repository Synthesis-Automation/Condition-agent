"""Condition-precedent association tests for web response enrichment."""

from __future__ import annotations

from app.web_api.condition_precedents import attach_condition_precedents


def test_condition_precedents_keep_reactions_with_their_own_evidence() -> None:
    payload = {
        "recommendations": [
            {
                "precedent_reaction_ids": ["reaction:1", "reaction:2"],
                "precedent_reaction_smiles": ["A>>B", "C>>D"],
                "precedent_reference_ids": ["REF1:second", "REF1:first"],
                "precedent_reaction_contexts": [
                    {
                        "reaction_id": "reaction:1",
                        "observation_id": "observation:1",
                        "reaction_smiles": "A>>B",
                        "reference_id": "REF1:first",
                    },
                    {
                        "reaction_id": "reaction:2",
                        "observation_id": "observation:2",
                        "reaction_smiles": "C>>D",
                        "reference_id": "REF1:second",
                    },
                ],
            }
        ]
    }
    references = {
        "REF1:first": {"reference_id": "REF1:first", "raw_reference": "First"},
        "REF1:second": {
            "reference_id": "REF1:second",
            "raw_reference": "Second",
        },
    }
    experiments = {
        "observation:observation:1": {
            "observation_id": "observation:1",
            "reaction_id": "reaction:1",
            "procedure_text": "Procedure one.",
        },
        "observation:observation:2": {
            "observation_id": "observation:2",
            "reaction_id": "reaction:2",
            "procedure_text": "Procedure two.",
        },
    }

    enriched = attach_condition_precedents(payload, references, experiments)
    precedents = enriched["recommendations"][0]["condition_precedents"]

    assert [value["reaction_smiles"] for value in precedents] == ["A>>B", "C>>D"]
    assert [
        value["reference_record"]["reference_id"] for value in precedents
    ] == ["REF1:first", "REF1:second"]
    assert [
        value["experimental_detail"]["procedure_text"] for value in precedents
    ] == ["Procedure one.", "Procedure two."]


def test_legacy_condition_precedents_do_not_guess_reference_pairings() -> None:
    payload = {
        "recommendations": [
            {
                "precedent_reaction_ids": ["reaction:1"],
                "precedent_reaction_smiles": ["A>>B"],
                "precedent_reference_ids": ["REF1:unpaired"],
            }
        ]
    }

    enriched = attach_condition_precedents(
        payload,
        {
            "REF1:unpaired": {
                "reference_id": "REF1:unpaired",
                "raw_reference": "Not row-associated",
            }
        },
        {},
    )

    precedent = enriched["recommendations"][0]["condition_precedents"][0]
    assert precedent["reaction_smiles"] == "A>>B"
    assert precedent["reference_id"] == ""
    assert precedent["reference_record"] is None
