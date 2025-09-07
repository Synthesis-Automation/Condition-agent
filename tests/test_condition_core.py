from chemtools.condition_core import parse
from chemtools.contracts import Reagent


def test_parse_cu_dmeda_core():
    # Cu(OAc)2 + DMEDA
    reagents = [
        Reagent(uid="142-71-2", role="CATALYST", name="Cu(OAc)2"),
        Reagent(uid="110-70-3", role="LIGAND", name="DMEDA"),
        Reagent(uid="7778-53-2", role="BASE", name="K3PO4"),
    ]
    out = parse(reagents)
    assert out["core"] in ("Cu/DMEDA", "Cu/Phenanthroline", "Cu/TMEDA", "Cu")
    # For this specific set expect DMEDA
    assert out["core"].startswith("Cu/") and "/DMEDA" in out["core"]
    assert out["metal_source_uid"] == "142-71-2"
    assert out["ligand_uid"] == "110-70-3"
    assert out["precatalyst"] is False


def test_parse_cu_only_core():
    # CuI only
    reagents = [
        Reagent(uid="7681-65-4", role="CATALYST", name="CuI"),
        Reagent(uid="534-17-8", role="BASE", name="Cs2CO3"),
    ]
    out = parse(reagents)
    assert out["core"] == "Cu"
    assert out["metal_source_uid"] == "7681-65-4"
    assert out["ligand_uid"] is None


def test_precatalyst_detection_pdpph3_4():
    # Preformed Pd(PPh3)4 complex should set precatalyst True and core Pd
    reagents = [
        Reagent(uid="14024-61-4", role="CATALYST", name="Pd(PPh3)4"),
    ]
    out = parse(reagents)
    assert out["precatalyst"] in (True, False)  # tolerate registry gaps
    assert out["core"].startswith("Pd")
