from chemtools.smiles import normalize_reaction

def test_normalize_reaction_smiles_basic():
    rsmi = "Brc1ccc(F)cc1.Nc1ccccc1>>"
    out = normalize_reaction(rsmi)
    assert 'normalized' in out and '>' in out['normalized']
    assert len(out['reactants']) == 2
    # Should not error on reasonable inputs
    assert out.get('errors') == []

