"""Regression tests for graph-derived departing-fragment identities."""

from dataclasses import asdict

import pytest
import reactive_taxonomy.departing_fragments as fragment_module
import reactive_taxonomy.reaction_parser as parser_module

from reactive_taxonomy import departing_fragment_tokens, featurize_reaction


def _tokens(reaction_smiles: str) -> tuple[str, ...]:
    analysis = featurize_reaction(reaction_smiles)
    assert analysis.reaction_signature is not None
    return departing_fragment_tokens(
        reaction_smiles,
        asdict(analysis.reaction_signature),
    )


def test_departing_fragment_tokens_separate_activation_states() -> None:
    alcohol = _tokens(
        "CCCCCCC(C)O.O=C1NC(=O)c2ccccc21.[K]>>CCCCCCC(C)N1C(=O)c2ccccc2C1=O"
    )
    mesylate = _tokens(
        "CCCCCCC(C)OS(C)(=O)=O.O=C1NC(=O)c2ccccc21.[K]>>CCCCCCC(C)N1C(=O)c2ccccc2C1=O"
    )
    tosylate = _tokens(
        "CCCCCCC(C)OS(=O)(=O)c1ccccc1.O=C1NC(=O)c2ccccc21.[K]"
        ">>CCCCCCC(C)N1C(=O)c2ccccc2C1=O"
    )
    substituted_tosylate = _tokens(
        "CCCCCCC(C)OS(=O)(=O)c1ccc(C)cc1.O=C1NC(=O)c2ccccc21.[K]"
        ">>CCCCCCC(C)N1C(=O)c2ccccc2C1=O"
    )

    assert alcohol and mesylate and tosylate
    assert len({alcohol, mesylate, tosylate}) == 3
    assert substituted_tosylate == tosylate


@pytest.mark.parametrize("leaving", ["O", "OS(C)(=O)=O", "OS(=O)(=O)c1ccccc1"])
def test_fragment_identity_matches_annotated_parsing_without_using_annotations(
    leaving: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    reaction = (
        f"CCCCCCC(C){leaving}.O=C1NC(=O)c2ccccc21.[K]>>CCCCCCC(C)N1C(=O)c2ccccc2C1=O"
    )
    signature = asdict(featurize_reaction(reaction).reaction_signature)
    fast = departing_fragment_tokens(reaction, signature)
    original = parser_module.parse_reaction_smiles

    def annotated(value: str, **kwargs):
        return original(value, include_molecular_interpretation=True)

    with monkeypatch.context() as patch:
        patch.setattr(fragment_module, "parse_reaction_smiles", annotated)
        assert departing_fragment_tokens(reaction, signature) == fast

    def unexpected_annotation(*args, **kwargs):
        raise AssertionError(
            "Graph-only fragment extraction must not annotate molecules"
        )

    monkeypatch.setattr(
        parser_module, "interpret_parsed_molecules", unexpected_annotation
    )
    assert departing_fragment_tokens(reaction, signature) == fast
    assert fast
