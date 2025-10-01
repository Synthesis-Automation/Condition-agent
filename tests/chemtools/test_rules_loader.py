from __future__ import annotations

from chemtools.rules import api


def test_load_crl_resolves_buchwald_reference() -> None:
    crl = api.load_default_crl()

    assert "Buchwald_CN" in crl.get("families", {})
    buchwald = crl["families"]["Buchwald_CN"]

    assert "defaults" in buchwald
    assert buchwald["defaults"]["pd_source"] == "Pd2(dba)3"

    playbook_ids = {pb.get("id") for pb in buchwald.get("playbooks", [])}
    assert "BH-ARCL-PRIM-ANILINE-GENERAL" in playbook_ids


def test_load_crl_includes_amide_formation() -> None:
    crl = api.load_default_crl()

    assert "Amide_Formation" in crl.get("families", {})
    amide = crl["families"]["Amide_Formation"]

    playbooks = amide.get("playbooks", [])
    assert playbooks, "Expected amide formation playbooks"

    ids = {pb.get("id") for pb in playbooks}
    assert "AMIDE-EDC-OXYMA-GENERAL" in ids
    assert "AM-BORON-B-C6F5-3" in ids


def test_load_crl_includes_suzuki_coupling() -> None:
    crl = api.load_default_crl()

    assert "Suzuki_Coupling" in crl.get("families", {})
    suzuki = crl["families"]["Suzuki_Coupling"]

    playbooks = suzuki.get("playbooks", [])
    assert len(playbooks) >= 10

    ids = {pb.get("id") for pb in playbooks}
    assert "SUZ-ARYL-I-BR-GENERAL-PPH3" in ids
    assert "SUZ-PROCESS-NI-MICELLAR" in ids


def test_load_crl_includes_ullmann_cn() -> None:
    crl = api.load_default_crl()

    assert "Ullmann_CN" in crl.get("families", {})
    ullmann = crl["families"]["Ullmann_CN"]

    defaults = ullmann.get("defaults", {})
    assert defaults.get("cu_source") == "CuI"

    playbook_ids = {pb.get("id") for pb in ullmann.get("playbooks", [])}
    assert "UL-ARBRI-PRIM-ANILINE-GENERAL" in playbook_ids

    guard_ids = {guard.get("id") for guard in ullmann.get("guards", [])}
    assert "UL-GUARD-ARCL-REQUIRES-2G-LIG" in guard_ids
