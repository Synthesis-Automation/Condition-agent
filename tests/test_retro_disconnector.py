from chemtools.retro.disconnector import RetronMatch, rank_disconnections


def test_rank_disconnections_prefers_concrete_precursors_over_conceptual_only(monkeypatch) -> None:
    monkeypatch.setattr(
        "chemtools.retro.disconnector.find_retrons",
        lambda mol, max_difficulty=0.8: [  # noqa: ARG005
            RetronMatch(
                retron_name="conceptual",
                reaction_name="easy_but_empty",
                difficulty=0.1,
                description="conceptual",
                notes="",
                precursor_hints=[],
                match_count=1,
                product_smarts="[*:1]-[*:2]",
            ),
            RetronMatch(
                retron_name="concrete",
                reaction_name="concrete_route",
                difficulty=0.4,
                description="concrete",
                notes="",
                precursor_hints=[],
                match_count=1,
                product_smarts="[*:1]-[*:2]",
            ),
        ],
    )
    monkeypatch.setattr(
        "chemtools.retro.retron_patterns.get_retron_by_name",
        lambda name: {"name": name},
    )

    def _fake_apply(mol, retron):  # noqa: ARG001
        if retron["name"] == "conceptual":
            return []
        return [("Brc1ccccc1", "Nc1ccccc1")]

    monkeypatch.setattr("chemtools.retro.disconnector._apply_retro_transforms", _fake_apply)

    out = rank_disconnections("c1ccccc1N", top_k=2)

    assert len(out) == 2
    assert out[0].reaction_name == "concrete_route"
    assert out[0].precursor_1 == "Brc1ccccc1"
    assert out[0].precursor_2 == "Nc1ccccc1"
    assert out[1].reaction_name == "easy_but_empty"
    assert out[1].precursor_1 == ""
    assert out[1].precursor_2 == ""
