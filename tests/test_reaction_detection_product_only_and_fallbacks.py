from __future__ import annotations

from chemtools.featurizers.formatters.detection_validation import (
    validate_detection_with_crk_key,
)
from chemtools.featurizers.formatters.reaction import _infer_multi_event_fallback_label


def test_product_only_family_can_match_when_reactants_missing() -> None:
    # CRK with no reacted motifs but formed vinylic halide should map to
    # Halogenation_unsaturated via product-only taxonomy slot matching.
    result = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|[] -> Alkenyl-Br",
        include_evidence=False,
    )
    assert result.get("reaction_type") == "Halogenation_unsaturated"


def test_event_fallback_label_from_single_event_kind() -> None:
    label = _infer_multi_event_fallback_label(
        {
            "events": [{"kind": "c_c_bond_formation", "confidence": 0.9}],
            "reaction_key_quality": {"reasons": []},
        }
    )
    assert label == "Event:C-C"


def test_event_fallback_label_for_no_motif_evidence_bucket() -> None:
    label = _infer_multi_event_fallback_label(
        {
            "events": [],
            "reaction_key_quality": {
                "reasons": [
                    "missing_formed_bond_and_product_motif_evidence",
                    "missing_bond_key",
                ]
            },
        }
    )
    assert label == "Unresolved:NoMotifEvidence"
