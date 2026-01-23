from chemtools.featurizers.spectator_rank import rank_spectator_groups


def test_rank_spectator_groups_basic_order() -> None:
    groups = ["OH", "Pyridine", "SH", "NH2", "CO2H", "F"]
    ranked = rank_spectator_groups(groups)
    assert ranked == ["Pyridine", "NH2", "CO2H", "OH", "SH", "F"]
