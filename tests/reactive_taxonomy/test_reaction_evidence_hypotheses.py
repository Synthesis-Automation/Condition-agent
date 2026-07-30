"""Evidence-provider and ambiguous-edit contracts."""

from reactive_taxonomy import featurize_reaction


FISCHER_REACTION = (
    "O=C1CCCCC1.Cl.NNc1ccc(F)cc1"
    ">>Fc1ccc2[nH]c3c(c2c1)CCCC3"
)


def test_ambiguous_correspondence_retains_typed_edit_hypotheses() -> None:
    analysis = featurize_reaction(FISCHER_REACTION)

    assert analysis.evidence_quality == "ambiguous_atom_correspondence"
    assert analysis.reaction_signature is None
    assert len(analysis.edit_hypotheses) == 2
    assert {
        candidate.provider for candidate in analysis.evidence_candidates
    } == {
        "exact_reconstruction",
        "exact_multi_event_reconstruction",
        "global_correspondence",
    }
    correspondence = next(
        candidate
        for candidate in analysis.evidence_candidates
        if candidate.provider == "global_correspondence"
    )
    assert correspondence.status == "ambiguous"
    assert correspondence.edit_hypotheses == analysis.edit_hypotheses

    for hypothesis in analysis.edit_hypotheses:
        assert hypothesis.hypothesis_id.startswith("REH1:")
        assert hypothesis.provider == "global_correspondence"
        assert hypothesis.evidence == "global_atom_correspondence"
        assert hypothesis.confidence == 0.8
        assert hypothesis.correspondence_count == 4
        assert hypothesis.edit_cost == (10, 6, 3)
        assert len(hypothesis.edits) == 9
        assert hypothesis.topology is not None
        assert hypothesis.topology.ring_count_delta == 1
        assert hypothesis.warnings == ("UNVERIFIED_EDIT_HYPOTHESIS",)


def test_edit_hypothesis_ids_are_deterministic_and_partner_order_invariant() -> None:
    product = "Fc1ccc2[nH]c3c(c2c1)CCCC3"
    forward = featurize_reaction(FISCHER_REACTION)
    repeated = featurize_reaction(FISCHER_REACTION)
    reversed_partners = featurize_reaction(
        f"NNc1ccc(F)cc1.Cl.O=C1CCCCC1>>{product}"
    )

    expected = tuple(
        hypothesis.hypothesis_id
        for hypothesis in forward.edit_hypotheses
    )
    assert expected
    assert expected == tuple(
        hypothesis.hypothesis_id
        for hypothesis in repeated.edit_hypotheses
    )
    assert expected == tuple(
        hypothesis.hypothesis_id
        for hypothesis in reversed_partners.edit_hypotheses
    )


def test_verified_reaction_exposes_provider_evidence_without_hypotheses() -> None:
    analysis = featurize_reaction("CCBr.N>>CCN")

    assert analysis.reaction_signature is not None
    assert not analysis.edit_hypotheses
    assert any(
        candidate.provider == "exact_reconstruction"
        and candidate.status == "verified"
        for candidate in analysis.evidence_candidates
    )
