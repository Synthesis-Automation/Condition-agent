"""Focused generation precedes validation limits and preserves mapped site choices."""

from dataclasses import asdict

import pytest

from core_retrosynthesis import build_generic_library, disconnect_generic_target_detailed
from core_retrosynthesis.generic_search import _apply, disconnect_operator_ladder_detailed
from tests.core_retrosynthesis_tests.test_generic_diverse_retrosynthesis import _row_from_reaction


@pytest.fixture(scope="module")
def library():
    reactions = [
        "[CH3:1][CH2:2]Br.[NH2:3][CH3:4]>>[CH3:1][CH2:2][NH:3][CH3:4]",
        "[CH3:4]Br.[CH3:1][CH2:2][NH2:3]>>[CH3:1][CH2:2][NH:3][CH3:4]",
        "[CH3:1][CH2:2]Br.[NH2:3][CH2:4][CH3:5]>>[CH3:1][CH2:2][NH:3][CH2:4][CH3:5]",
    ]
    return build_generic_library(tuple(
        _row_from_reaction(reaction, reaction_id=f"r-{i}", reference_id=f"p-{i}")
        for i, reaction in enumerate(reactions)
    ), engine="reaction_core", admission_mode="data_driven", levels=("L0", "L1", "L2"))


@pytest.mark.parametrize("bond,expected", [([1, 2], "CCBr.CN"), ([2, 3], "CBr.CCN")])
def test_two_sites_can_each_be_selected_under_one_validation_budget(library, bond, expected):
    candidates, diagnostics = disconnect_generic_target_detailed(
        "CCNC", library, top_k=1, max_candidates_to_validate=1,
        required_disconnection_bond=bond, focus_target_smiles="CCNC",
    )
    assert candidates, diagnostics
    assert candidates[0].precursor_smiles == expected
    assert candidates[0].bond_focus_check.status == "verified"
    assert diagnostics.validation_attempt_count == 1
    assert diagnostics.focus_rejected_count > 0


def test_reordered_target_and_unrestricted_default_are_deterministic(library):
    normal, counts = disconnect_generic_target_detailed("CCNC", library)
    explicit, explicit_counts = disconnect_generic_target_detailed(
        "CCNC", library, required_disconnection_bond=None, focus_target_smiles=None,
    )
    assert normal == explicit and counts == explicit_counts
    options = {"required_disconnection_bond": [1, 2], "focus_target_smiles": "CCNC"}
    left, _ = disconnect_generic_target_detailed("CCNC", library, **options)
    right, _ = disconnect_generic_target_detailed("CNCC", library, **options)
    assert left == right


def test_symmetric_sites_survive_precursor_deduplication(library):
    for bond in ([1, 2], [2, 3]):
        values, _ = disconnect_generic_target_detailed(
            "CCNCC", library, required_disconnection_bond=bond, focus_target_smiles="CCNCC",
        )
        assert values
        assert all(value.bond_focus_check.target_atom_ids == tuple(bond) for value in values)
        assert "CCBr.CCN" in {value.precursor_smiles for value in values}
    template = next(item for item in library.templates if item.abstraction_level == "L0")
    outcomes = _apply(template.reaction_smarts, "CCNCC", preserve_correspondences=True)
    assert len(outcomes) > len({precursors for precursors, _ in outcomes})


def test_no_qualifying_bond_is_not_silently_relaxed(library):
    values, diagnostics = disconnect_operator_ladder_detailed(
        "CCNC", library, required_disconnection_bond=[0, 1], focus_target_smiles="CCNC",
    )
    assert values == ()
    assert diagnostics.validation_attempt_count == 0
    assert sum(item.focus_rejected_count for _, item in diagnostics.level_diagnostics) > 0
    assert asdict(diagnostics.level_diagnostics[0][1])["focus_rejection_examples"]


def test_failed_final_focus_confirmation_cannot_be_returned(library, monkeypatch):
    import core_retrosynthesis.generic_search as search
    from dataclasses import replace

    real = search.check_disconnection_bond

    def conflicting(*args, **kwargs):
        check = real(*args, **kwargs)
        return replace(check, status="unresolved") if "observation" in kwargs else check

    monkeypatch.setattr(search, "check_disconnection_bond", conflicting)
    values, diagnostics = disconnect_generic_target_detailed(
        "CCNC", library, required_disconnection_bond=[1, 2], focus_target_smiles="CCNC",
    )
    assert not values and diagnostics.focus_validation_rejected_count > 0


def test_controlled_effectiveness_at_equal_template_and_validation_budgets(library):
    baseline, baseline_counts = disconnect_generic_target_detailed(
        "CCNC", library, top_k=1, max_templates_to_apply=40, max_candidates_to_validate=1,
    )
    expected = {(1, 2): "CCBr.CN", (2, 3): "CBr.CCN"}
    baseline_hits = sum(any(c.precursor_smiles == value for c in baseline)
                        for value in expected.values())
    focused_hits = 0
    for bond, precursor in expected.items():
        values, counts = disconnect_generic_target_detailed(
            "CCNC", library, top_k=1, max_templates_to_apply=40, max_candidates_to_validate=1,
            required_disconnection_bond=bond, focus_target_smiles="CCNC",
        )
        focused_hits += any(c.precursor_smiles == precursor for c in values)
        assert counts.applied_template_count == baseline_counts.applied_template_count
        assert counts.validation_attempt_count == baseline_counts.validation_attempt_count == 1
    assert baseline_hits == 1 and focused_hits == 2
