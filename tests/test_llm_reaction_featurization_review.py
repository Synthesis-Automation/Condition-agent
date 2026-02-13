from __future__ import annotations

from dataclasses import dataclass

from llmtools.reaction_featurization_review import (
    LLMReactionFeaturizationOptions,
    build_reaction_featurization_prompt,
    review_reaction_featurization,
)


@dataclass
class _FakeLLMResponse:
    content: str
    model: str = "gpt-test"
    prompt_tokens: int = 10
    completion_tokens: int = 20
    total_tokens: int = 30
    latency_ms: float = 12.5


class _FakeLLMClient:
    def __init__(self, content: str) -> None:
        self._content = content

    def chat(self, prompt: str, temperature: float = 0.0, max_tokens: int = 0):
        assert "Reaction SMILES" in prompt
        return _FakeLLMResponse(content=self._content)


def test_build_reaction_featurization_prompt_includes_context() -> None:
    prompt = build_reaction_featurization_prompt(
        {
            "reaction_smiles": "CCO.CN>>CCN",
            "deterministic_reaction_type": "Unknown",
            "deterministic_confidence": 0.2,
            "reacted_motifs": ["Ar-Br"],
            "formed_motifs": ["Ar-NR2"],
        }
    )
    assert "CCO.CN>>CCN" in prompt
    assert "Unknown" in prompt
    assert "Ar-Br" in prompt
    assert "Ar-NR2" in prompt
    assert "Stoichiometry delta" in prompt
    assert "Reaction event kinds" in prompt


def test_review_reaction_featurization_parses_json_fences() -> None:
    options = LLMReactionFeaturizationOptions(
        enabled=True,
        provider="openai",
        model="gpt-test",
    )
    fake = _FakeLLMClient(
        """```json
{"suggested_reaction_type":"C_N_Coupling","confidence":0.82,"rationale":"motif pattern fit","requires_human_review":false,"uncertainty_flags":["low_mapping_confidence"]}
```"""
    )

    out = review_reaction_featurization(
        {"reaction_smiles": "CCO.CN>>CCN"},
        options,
        client=fake,  # type: ignore[arg-type]
    )

    assert out["status"] == "ok"
    analysis = out.get("analysis") or {}
    assert analysis.get("suggested_reaction_type") == "C_N_Coupling"
    assert analysis.get("confidence") == 0.82
    assert analysis.get("requires_human_review") is False


def test_review_reaction_featurization_parses_extended_optional_fields() -> None:
    options = LLMReactionFeaturizationOptions(
        enabled=True,
        provider="openai",
        model="gpt-test",
    )
    fake = _FakeLLMClient(
        """{
  "suggested_reaction_type":"SNAr_CN",
  "confidence":0.73,
  "rationale":"C-Cl displacement with C-N formation",
  "requires_human_review":true,
  "uncertainty_flags":["taxonomy_boundary_case"],
  "mechanistic_family":"SNAr",
  "mechanistic_rationale":"Electron-poor heteroaryl C-Cl plus N nucleophile signal",
  "tautomer_or_representation_issue":true,
  "taxonomy_gap_suspected":true,
  "taxonomy_gap_note":"Hydrazine-like motif in product set not explicit in slot constraints",
  "deterministic_checks_used":["reaction_key","event_kinds","stoichiometry_delta"]
}"""
    )

    out = review_reaction_featurization(
        {"reaction_smiles": "CCO.CN>>CCN"},
        options,
        client=fake,  # type: ignore[arg-type]
    )

    assert out["status"] == "ok"
    analysis = out.get("analysis") or {}
    assert analysis.get("mechanistic_family") == "SNAr"
    assert analysis.get("tautomer_or_representation_issue") is True
    assert analysis.get("taxonomy_gap_suspected") is True
    assert analysis.get("deterministic_checks_used") == [
        "reaction_key",
        "event_kinds",
        "stoichiometry_delta",
    ]
