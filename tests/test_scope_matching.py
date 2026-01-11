from chemtools.HTE import recommender as r


def test_scope_expands_any_conh2_to_alkyl() -> None:
    motif_sets = r._load_motif_sets()
    scope_map = r._load_scope_map()
    expanded = r._expand_motif_tokens(["Any-CONH2"], motif_sets, scope_map)
    assert "Alkyl-CONH2" in expanded
