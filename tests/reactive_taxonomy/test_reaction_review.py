"""Regression tests for the shared reaction-review presentation contract."""

from reactive_taxonomy import (
    REACTION_REVIEW_SUMMARY_SCHEMA_VERSION,
    build_reaction_review_summary,
    featurize_reaction,
    format_reaction_review_summary,
)


def test_review_summary_combines_detailed_and_graphic_labels() -> None:
    analysis = featurize_reaction(
        "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]"
    )

    summary = build_reaction_review_summary(analysis)

    assert summary.schema_version == REACTION_REVIEW_SUMMARY_SCHEMA_VERSION
    assert "H₂C=CH₂ → H₃C–CH₃" in summary.detailed_reaction_label
    assert summary.reaction_core_equation == "C=C → C–C"
    assert "unsaturated partner" in summary.electronic_steric_analysis
    assert "access open" in summary.electronic_steric_analysis


def test_live_analysis_and_serialized_record_share_review_rendering() -> None:
    analysis = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    live = build_reaction_review_summary(analysis)
    serialized = analysis.to_dict()
    record = build_reaction_review_summary(serialized)

    assert record == live
    rendered = format_reaction_review_summary(live)
    assert rendered.splitlines()[0].startswith("Detailed reaction label:")
    assert rendered.splitlines()[1] == "Graphic core label: Unavailable"
    assert rendered.splitlines()[2] == "Spectators: None detected"
    assert rendered.splitlines()[3].startswith(
        "Electronic / steric analysis: electrophile:"
    )


def test_review_summary_formats_duplicate_spectators_without_distances() -> None:
    summary = build_reaction_review_summary(
        {
            "reaction_signature": {
                "spectator_groups": [
                    {"group_id": "ether", "chemist_label": "R–O–R"},
                    {"group_id": "ether", "chemist_label": "R–O–R"},
                ],
                "partners": [],
            }
        }
    )

    assert summary.spectators == "2× R–O–R [ether]"


def test_review_summary_shows_general_motif_and_specific_limiter() -> None:
    analysis = featurize_reaction(
        "O=[C:15](O)[CH2:2][CH3:1]."
        "[cH:3]1[cH:4][c:5]([CH3:6])[n:7][c:8]2[cH:9][cH:10]"
        "[c:11]([F:12])[cH:13][c:14]12"
        ">>[CH3:1][CH2:2][c:3]1[cH:4][c:5]([CH3:6])[n:7]"
        "[c:8]2[cH:9][cH:10][c:11]([F:12])[cH:13][c:14]12"
    )

    summary = build_reaction_review_summary(analysis)
    rendered = format_reaction_review_summary(summary)

    assert summary.reaction_core_equation == (
        "R–C(=O)OH + Ar–H → R–Ar"
    )
    assert summary.core_limiter == (
        "R = R′–CH₂– (primary alkyl); Ar = HetAr"
    )
    assert summary.atom_level_core_equation == (
        "C(H)2(R)(C(=O)(O-H)) + C(H)(:HetAr)2 → C(H)2(R)(ArC)"
    )
    assert "Graphic core limiter:" in rendered
    assert "Atom-level core:" in rendered
