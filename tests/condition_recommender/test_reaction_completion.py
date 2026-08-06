from dataclasses import replace

import pytest

from condition_recommender import (
    build_completion_selection,
    build_completed_reaction_smiles,
    propose_reaction_completion,
)


_COUNTERION_AZIDATION_QUERY = (
    "C=Cn1cc[n+](Cc2c(F)c(F)c(F)c(F)c2F)c1.[Br-]"
    ">>C=Cn1cc[n+](Cc2c(F)c(F)c(N=[N+]=[N-])c(F)c2F)c1.[Br-]"
)


def test_incomplete_azidation_proposes_curated_source_choices() -> None:
    proposal = propose_reaction_completion(_COUNTERION_AZIDATION_QUERY)

    assert proposal.status == "confirmation_recommended"
    assert proposal.provenance == "system_proposed"
    assert proposal.query_reaction_smiles == _COUNTERION_AZIDATION_QUERY
    assert len(proposal.requirements) == 1
    requirement = proposal.requirements[0]
    assert requirement.rooted_fragment_smiles == "*N=[N+]=[N-]"
    assert [value.option_kind for value in requirement.options] == [
        "compatible_source_class",
        "registered_substance",
        "registered_substance",
        "unresolved",
    ]
    assert {
        value.substance_id
        for value in requirement.options
        if value.substance_id
    } == {"cas:26628-22-8", "cas:4648-54-8"}


def test_edited_known_identifier_resolves_without_rewriting_query() -> None:
    proposal = propose_reaction_completion(_COUNTERION_AZIDATION_QUERY)
    requirement = proposal.requirements[0]

    selection = build_completion_selection(
        proposal,
        requirement.requirement_id,
        custom_identifier="26628-22-8",
    )

    assert selection.selection_kind == "registered_substance"
    assert selection.provenance == "user_edited"
    assert selection.substance_id == "cas:26628-22-8"
    assert selection.raw_identifier == "26628-22-8"
    assert proposal.query_reaction_smiles == _COUNTERION_AZIDATION_QUERY


def test_confirmed_sodium_azide_builds_separate_completed_query() -> None:
    proposal = propose_reaction_completion(_COUNTERION_AZIDATION_QUERY)
    requirement = proposal.requirements[0]
    sodium_azide = next(
        value
        for value in requirement.options
        if value.substance_id == "cas:26628-22-8"
    )
    selection = build_completion_selection(
        proposal,
        requirement.requirement_id,
        option_id=sodium_azide.option_id,
    )

    completed, warnings = build_completed_reaction_smiles(
        _COUNTERION_AZIDATION_QUERY,
        (selection,),
    )

    assert completed is not None
    assert ">[Na+].[N-]=[N+]=[N-]>" in completed
    assert completed != _COUNTERION_AZIDATION_QUERY
    assert warnings == ()
    assert proposal.query_reaction_smiles == _COUNTERION_AZIDATION_QUERY


def test_unknown_edited_source_remains_explicitly_unresolved() -> None:
    proposal = propose_reaction_completion(_COUNTERION_AZIDATION_QUERY)
    requirement = proposal.requirements[0]

    selection = build_completion_selection(
        proposal,
        requirement.requirement_id,
        custom_identifier="unregistered azide source",
    )

    assert selection.selection_kind == "custom_identifier"
    assert selection.provenance == "user_edited"
    assert not selection.resolved
    assert selection.substance_id is None


def test_stale_completion_selection_is_rejected() -> None:
    proposal = propose_reaction_completion(_COUNTERION_AZIDATION_QUERY)
    requirement = proposal.requirements[0]
    option = requirement.options[0]
    selection = build_completion_selection(
        proposal,
        requirement.requirement_id,
        option_id=option.option_id,
    )

    with pytest.raises(ValueError, match="different proposal"):
        from condition_recommender.reaction_completion import (
            validate_completion_selections,
        )

        validate_completion_selections(
            proposal,
            (replace(selection, proposal_id="RCP1:stale"),),
        )
