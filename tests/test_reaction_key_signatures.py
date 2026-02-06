from chemtools.recommend.recommender import (
    _reaction_key_to_signatures,
    _reactant_types_to_signature,
    _signature_matches,
)


def test_reaction_key_to_signatures_core_ext():
    key = (
        "|Ar-Br|R-NH2 -> Ar-NR2 | bond_formed: C(ar)-N | bond_broken: Br-C(ar) | spectators: Pyridine"
    )
    core, ext = _reaction_key_to_signatures(key)
    assert core == "Ar-Br|R-NH2"
    assert ext == "Ar-Br|Pyridine|R-NH2"


def test_reaction_key_to_signatures_empty():
    core, ext = _reaction_key_to_signatures("")
    assert core == ""
    assert ext == ""


def test_reactant_types_to_signature():
    sig = _reactant_types_to_signature(["R-NH2", "Ar-Br"])
    assert sig == "Ar-Br|R-NH2"


def test_signature_matches_aromatic_scaffold_equivalents():
    query_sets = [set(["Ar-B(OH)2", "HeteroAr-Br"])]
    assert _signature_matches(query_sets, "Ar-B(OH)2|Ar-Br") is True
