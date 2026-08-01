"""Review rendering uses the canonical reaction label and minimum core."""

from reactive_taxonomy import featurize_reaction
from reactive_taxonomy.reaction_review import build_reaction_review_summary


def test_review_summary_uses_single_label_contract() -> None:
    analysis = featurize_reaction(
        "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    )
    summary = build_reaction_review_summary(analysis)
    assert analysis.reaction_label is not None
    assert summary.detailed_reaction_label == analysis.reaction_label.detailed
    assert analysis.reaction_core is not None
    assert summary.reaction_core_equation == analysis.reaction_core.generic_label


def test_unavailable_review_label_remains_explicit() -> None:
    summary = build_reaction_review_summary(featurize_reaction("CC.O>>N#N"))
    assert summary.detailed_reaction_label == (
        "No normalized edits or supported structural interpretation"
    )
