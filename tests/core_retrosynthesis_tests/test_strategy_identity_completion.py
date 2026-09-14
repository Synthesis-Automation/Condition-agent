"""Complete missing strategy identities without replacing existing identities."""

from core_retrosynthesis.generic_compiler import analyze_generic_reaction


def test_formed_bond_site_keeps_existing_namespace():
    identity = analyze_generic_reaction("Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
    assert identity is not None
    assert identity.disconnection_site_key.startswith("SITE1:")
    assert identity.synthon_signature.startswith("SYN1:")


def test_departure_only_site_and_aromatic_synthon_receive_graph_identities():
    reaction = "[N:1]#[C:2][c:3]1[n:4][n:5]([C:901]([CH3:900])=[O:902])[c:6]2[cH:7][cH:8][cH:9][cH:10][c:11]12>>[N:1]#[C:2][c:3]1[n:4][nH:5][c:6]2[cH:7][cH:8][cH:9][cH:10][c:11]12"
    identity = analyze_generic_reaction(reaction)
    assert identity is not None
    assert identity.disconnection_site_key.startswith("SITE2:")
    assert identity.synthon_signature.startswith("SYN2:")
