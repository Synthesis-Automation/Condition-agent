from chem_coworker.agent import ChemCoworker
from chem_coworker.tools import REGISTRY


def _make_agent() -> ChemCoworker:
    agent = object.__new__(ChemCoworker)
    agent.registry = REGISTRY
    return agent


def test_output_gate_redacts_invalid_reaction_smiles(monkeypatch) -> None:
    agent = _make_agent()
    monkeypatch.setattr(
        "chem_coworker.tools._helpers._validate_reaction_smiles",
        lambda reaction_smiles, require_product=True: (reaction_smiles, "invalid reactant"),  # noqa: ARG005
    )

    answer, warnings, penalty = agent._apply_output_verification_gate(
        answer="Reaction SMILES: `A.B>>C`",
        tool_results={},
        task_type="retrosynthesis",
    )

    assert "[INVALID_REACTION_SMILES]" in answer
    assert "output verification gate" in answer.lower()
    assert any("removed invalid reaction smiles" in w.lower() for w in warnings)
    assert penalty > 0.0


def test_output_gate_flags_unbacked_condition_claims() -> None:
    agent = _make_agent()
    answer, warnings, penalty = agent._apply_output_verification_gate(
        answer=(
            "Recommended conditions: catalyst Pd(OAc)2, ligand XPhos, "
            "base Cs2CO3, solvent dioxane, temperature 90 C."
        ),
        tool_results={},
        task_type="forward_synthesis",
    )

    assert any("condition recommendations" in w.lower() for w in warnings)
    assert "output verification gate" in answer.lower()
    assert penalty > 0.0


def test_output_gate_flags_unproven_route_steps(monkeypatch) -> None:
    agent = _make_agent()
    monkeypatch.setattr(
        "chem_coworker.tools._helpers._validate_reaction_smiles",
        lambda reaction_smiles, require_product=True: (reaction_smiles, None),  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.composite._evaluate_synthesis_proposal",
        lambda **kwargs: {  # noqa: ARG005
            "success": True,
            "verdict": "PASS",
            "critical_failures": [],
            "warnings": [],
        },
    )

    answer, warnings, penalty = agent._apply_output_verification_gate(
        answer="Step 1: A.B>>C\nStep 2: C.D>>E",
        tool_results={"plan_route": {"success": True, "route": []}},
        task_type="retrosynthesis",
    )

    assert any("multiple explicit reaction steps" in w.lower() for w in warnings)
    assert "route details include unverified explicit reaction steps" in answer.lower()
    assert penalty >= 0.1


def test_output_gate_keeps_supported_answer(monkeypatch) -> None:
    agent = _make_agent()
    monkeypatch.setattr(
        "chem_coworker.tools._helpers._validate_reaction_smiles",
        lambda reaction_smiles, require_product=True: (reaction_smiles, None),  # noqa: ARG005
    )
    monkeypatch.setattr(
        "chem_coworker.tools.composite._evaluate_synthesis_proposal",
        lambda **kwargs: {  # noqa: ARG005
            "success": True,
            "verdict": "PASS",
            "critical_failures": [],
            "warnings": [],
        },
    )

    answer, warnings, penalty = agent._apply_output_verification_gate(
        answer=(
            "Reaction SMILES: A.B>>C\n"
            "Recommended conditions: catalyst Pd(OAc)2, ligand XPhos, base Cs2CO3, solvent dioxane."
        ),
        tool_results={
            "recommend_reaction_conditions": {
                "success": True,
                "recommendations": [
                    {
                        "catalyst": "Pd(OAc)2",
                        "ligand": "XPhos",
                        "base": "Cs2CO3",
                        "solvent": "dioxane",
                        "temperature": 90,
                    }
                ],
            },
            "retrosynthesis_step": {
                "success": True,
                "top_disconnection": {"precursor_1": "A", "precursor_2": "B"},
            },
        },
        task_type="retrosynthesis",
    )

    assert warnings == []
    assert penalty == 0.0
    assert "output verification gate" not in answer.lower()
