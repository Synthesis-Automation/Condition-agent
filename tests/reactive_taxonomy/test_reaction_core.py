"""Minimum-core regressions for mapped and inferred structural evidence."""

from reactive_taxonomy import featurize_reaction


def test_mapped_reaction_builds_minimum_core() -> None:
    result = featurize_reaction("[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]")
    core = result.reaction_core
    assert core is not None
    assert core.evidence_status == "verified"
    assert core.generic_label
    assert core.shape_core_key.startswith("RSH2:")
    assert core.center_transition_key.startswith("RCS2:")


def test_inferred_correspondence_builds_core_with_internal_atom_ids() -> None:
    result = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    core = result.reaction_core
    assert core is not None
    assert core.evidence_status == "inferred"
    assert core.generic_label == "Ar–B + Ar–Br → Ar–Ar"
    assert "site:" not in " ".join(core.participant_tokens)


def test_small_split_reagent_is_handled_without_reaction_rule() -> None:
    result = featurize_reaction("C=C.BrBr>>BrCCBr")
    assert result.evidence_quality == "global_atom_correspondence"
    assert result.reaction_core is not None
    assert "Br–Br" in result.reaction_core.generic_label
    assert result.interpretation is not None
    assert result.interpretation.primary_pattern_id == "net_addition"


def test_chemically_distinct_atom_origin_hypotheses_do_not_force_one_core() -> None:
    result = featurize_reaction("CC(O)C.CO>>COC(C)C")
    assert result.evidence_quality == "ambiguous_atom_correspondence"
    assert len(result.edit_hypotheses) == 2
    assert result.reaction_core is None


def test_core_is_deterministic_and_partner_order_invariant() -> None:
    forward = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    reverse = featurize_reaction(
        "OB(O)c1ccccc1.Brc1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    assert forward.reaction_core is not None
    assert reverse.reaction_core is not None
    assert forward.reaction_core.shape_core_key == reverse.reaction_core.shape_core_key
    assert forward.reaction_core.generic_label == reverse.reaction_core.generic_label
