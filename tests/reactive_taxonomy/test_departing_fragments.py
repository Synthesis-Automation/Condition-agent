"""Regression tests for graph-derived departing-fragment identities."""

from dataclasses import asdict

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
        "CCCCCCC(C)O.O=C1NC(=O)c2ccccc21.[K]"
        ">>CCCCCCC(C)N1C(=O)c2ccccc2C1=O"
    )
    mesylate = _tokens(
        "CCCCCCC(C)OS(C)(=O)=O.O=C1NC(=O)c2ccccc21.[K]"
        ">>CCCCCCC(C)N1C(=O)c2ccccc2C1=O"
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
