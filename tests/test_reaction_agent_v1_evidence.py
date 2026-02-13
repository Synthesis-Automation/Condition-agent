"""Tests for canonical reaction evidence model."""

from reaction_agent_v1.evidence import ReactionEvidence


def test_evidence_merge_diff_payload_populates_core_sections() -> None:
    evidence = ReactionEvidence(reaction_smiles="A>>B")
    payload = {
        "analysis": {
            "reaction_smiles": "A>>B",
            "principal_pair": {"reactant_smiles": "A", "product_smiles": "B"},
            "mcs_ratio": 0.5,
            "core_formula_delta": {"N": 1},
            "side_formula_delta": {"N": 1},
            "reacted_motifs": ["Ar-Cl"],
            "formed_motifs": ["Ar-NH2"],
            "reaction_key": "demo_key",
            "taxonomy_candidates": [{"reaction_type": "C_N_Coupling", "deterministic_score": 0.8}],
            "decision": {"reaction_type": "C_N_Coupling", "confidence": 0.8},
        }
    }
    evidence.merge_diff_payload(payload)

    assert evidence.has_diff()
    assert evidence.diff["principal_pair"]["reactant_smiles"] == "A"
    assert evidence.detection["reaction_key"] == "demo_key"
    assert evidence.taxonomy_candidates[0]["reaction_type"] == "C_N_Coupling"
    assert evidence.provisional_decision["reaction_type"] == "C_N_Coupling"
    assert evidence.validation == {}
