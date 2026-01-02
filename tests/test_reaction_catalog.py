from chemtools.taxonomy.reaction_catalog import load_reaction_catalog


def test_reaction_catalog_expands_motif_sets() -> None:
    definitions, _ = load_reaction_catalog()

    suzuki = definitions["suzuki_miyaura"]
    assert "Ar-OTf" in suzuki.reactants["electrophiles"].allowed
    assert "Ar-B(OH)2" in suzuki.reactants["nucleophiles"].allowed

    sonogashira = definitions["sonogashira"]
    assert "Ar-OMs" in sonogashira.reactants["electrophiles"].allowed
    assert "Ar-F" not in sonogashira.reactants["electrophiles"].allowed

    snar = definitions["snar_cn"]
    assert "Ar-NO2" in snar.reactants["activators"].allowed
    assert "Pyridine" in snar.reactants["activators"].allowed

    cn = definitions["c_n_cross_coupling"]
    assert "Ar-NH2" in cn.reactants["nucleophiles"].allowed
    assert "Ar-CONHR" in cn.reactants["nucleophiles"].allowed
    assert "Indole" in cn.reactants["nucleophiles"].allowed

    sn2 = definitions["sn2_substitution"]
    assert "R-OTs" in sn2.reactants["electrophiles"].allowed

    reductive = definitions["reductive_amination"]
    assert "Any-CHO" in reductive.reactants["electrophiles"].allowed

    rcm = definitions["ring_closing_metathesis"]
    assert "Any-CHCH2" in rcm.reactants["nucleophiles"].allowed
