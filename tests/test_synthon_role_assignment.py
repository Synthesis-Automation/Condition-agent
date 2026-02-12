from chemtools.precedent.loader import _pick_electrophile_nucleophile
from chemtools.recommend.utils import pick_electrophile_nucleophile
from chemtools.synthon import classify_reactant_synthons, select_electrophile_nucleophile


def test_classify_reactant_synthons_sp2_electrophile() -> None:
    hits = classify_reactant_synthons("Brc1ccccc1")
    assert any(hit.role == "electrophile" for hit in hits)
    assert any(hit.synthon_id == "sp2_electrophile" for hit in hits)


def test_classify_reactant_synthons_organoboron_partner() -> None:
    hits = classify_reactant_synthons("OB(O)c1ccccc1")
    assert any(hit.role == "nucleophile" for hit in hits)
    assert any(hit.synthon_id == "organoboron_partner" for hit in hits)


def test_select_pair_cn_coupling() -> None:
    elec, nuc = select_electrophile_nucleophile(["Nc1ccccc1", "Brc1ccccc1"])
    assert elec == "Brc1ccccc1"
    assert nuc == "Nc1ccccc1"


def test_select_pair_suzuki_reverse_input_order() -> None:
    elec, nuc = select_electrophile_nucleophile(["OB(O)c1ccccc1", "Brc1ccccc1"])
    assert elec == "Brc1ccccc1"
    assert nuc == "OB(O)c1ccccc1"


def test_select_pair_prefers_distinct_reactants_for_ambifunctional_substrate() -> None:
    elec, nuc = select_electrophile_nucleophile(["Nc1ccc(Br)cc1", "OB(O)c1ccccc1"])
    assert elec == "Nc1ccc(Br)cc1"
    assert nuc == "OB(O)c1ccccc1"


def test_select_pair_single_reactant_fallback() -> None:
    elec, nuc = select_electrophile_nucleophile(["Brc1ccccc1"])
    assert elec == "Brc1ccccc1"
    assert nuc == ""


def test_wrapper_functions_delegate_to_synthon_selector() -> None:
    reactants = ["OB(O)c1ccccc1", "Brc1ccccc1"]
    expected = select_electrophile_nucleophile(reactants)
    assert pick_electrophile_nucleophile(reactants) == expected
    assert _pick_electrophile_nucleophile(reactants) == expected
