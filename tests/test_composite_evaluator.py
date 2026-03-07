from chem_coworker.tools.composite import (
    _evaluate_synthesis_proposal,
    _validate_synthesis_proposal,
)
from chem_coworker.tools import REGISTRY


def test_evaluate_synthesis_proposal_reaction_mode_adds_score_and_grade(monkeypatch):
    monkeypatch.setattr(
        "chem_coworker.tools.reaction_eval._evaluate_reaction",
        lambda reaction_smiles, reaction_type="": {  # noqa: ARG005
            "success": True,
            "verdict": "PASS",
            "critical_failures": [],
            "warnings": [],
        },
    )
    monkeypatch.setattr(
        "chem_coworker.tools.chemistry._detect_reaction_type",
        lambda reaction_smiles: {  # noqa: ARG005
            "success": True,
            "reaction_type": "suzuki_miyaura",
            "reaction_type_id": "suzuki_miyaura",
            "confidence": 0.91,
        },
    )
    monkeypatch.setattr(
        "chemtools.taxonomy.reaction_catalog.resolve_reaction_type",
        lambda value: str(value or "").strip().lower(),
    )

    out = _evaluate_synthesis_proposal(
        mode="reaction",
        reaction_smiles="A.B>>C",
        reaction_name="suzuki_miyaura",
    )
    assert out["success"] is True
    assert out["evaluation_score"] >= 85
    assert out["evaluation_grade"] == "A"
    assert out["evaluation_summary"]["consistency"]["checked"] is True
    assert out["evaluation_summary"]["consistency"]["matches_expected_reaction"] is True


def test_evaluate_synthesis_proposal_reaction_mode_penalizes_type_mismatch(monkeypatch):
    monkeypatch.setattr(
        "chem_coworker.tools.reaction_eval._evaluate_reaction",
        lambda reaction_smiles, reaction_type="": {  # noqa: ARG005
            "success": True,
            "verdict": "PASS",
            "critical_failures": [],
            "warnings": [],
        },
    )
    monkeypatch.setattr(
        "chem_coworker.tools.chemistry._detect_reaction_type",
        lambda reaction_smiles: {  # noqa: ARG005
            "success": True,
            "reaction_type": "suzuki_miyaura",
            "reaction_type_id": "suzuki_miyaura",
            "confidence": 0.88,
        },
    )
    monkeypatch.setattr(
        "chemtools.taxonomy.reaction_catalog.resolve_reaction_type",
        lambda value: str(value or "").strip().lower(),
    )

    out = _evaluate_synthesis_proposal(
        mode="reaction",
        reaction_smiles="A.B>>C",
        reaction_name="buchwald_hartwig",
    )
    assert out["success"] is True
    assert out["evaluation_score"] <= 80
    assert out["evaluation_summary"]["consistency"]["matches_expected_reaction"] is False


def test_evaluate_synthesis_proposal_retro_penalizes_no_simplification(monkeypatch):
    monkeypatch.setattr(
        "chem_coworker.tools.reaction_eval._check_retro_consistency",
        lambda product_smiles, precursor_1, precursor_2, reaction_name="": {  # noqa: ARG005
            "success": True,
            "verdict": "PASS",
            "critical_failures": [],
            "warnings": [],
            "complexity_check": {"simplification_achieved": False},
        },
    )

    out = _evaluate_synthesis_proposal(
        mode="retro",
        product_smiles="P",
        precursor_1="A",
        precursor_2="B",
    )
    assert out["success"] is True
    assert out["evaluation_score"] <= 70
    assert out["evaluation_grade"] in {"B", "C", "D", "F"}


def test_validate_synthesis_proposal_backward_alias_now_returns_eval_fields(monkeypatch):
    monkeypatch.setattr(
        "chem_coworker.tools.reaction_eval._evaluate_reaction",
        lambda reaction_smiles, reaction_type="": {  # noqa: ARG005
            "success": True,
            "verdict": "PASS_WITH_WARNINGS",
            "critical_failures": [],
            "warnings": ["warn"],
        },
    )
    monkeypatch.setattr(
        "chem_coworker.tools.chemistry._detect_reaction_type",
        lambda reaction_smiles: {"success": False},
    )

    out = _validate_synthesis_proposal(mode="reaction", reaction_smiles="A>>B")
    assert out["success"] is True
    assert "evaluation_score" in out
    assert "evaluation_grade" in out


def test_registry_contains_dedicated_evaluator_tool():
    plugin = REGISTRY._plugins.get("evaluate_synthesis_proposal")
    assert plugin is not None
    assert plugin.llm_exposed is True
