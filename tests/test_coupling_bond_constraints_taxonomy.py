from __future__ import annotations

from chemtools.reaction.feasibility import (
    validate_detection_with_crk_key,
)


def test_c_o_coupling_requires_c_o_bond_formation_in_taxonomy_constraints() -> None:
    ok = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|Ar-Cl|Ar-OH -> Ar-OR | bond_formed: C-O | bond_broken: C-Cl",
    )
    bad = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|Ar-Cl|Ar-OH -> Ar-OR | bond_formed: C-C | bond_broken: C-Cl",
    )
    assert ok["reaction_type"] == "C_O_Coupling"
    assert bad["reaction_type"] == "Unknown"


def test_c_s_coupling_requires_c_s_bond_formation_in_taxonomy_constraints() -> None:
    ok = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|Ar-Cl|Ar-SH -> Ar-SR | bond_formed: C-S | bond_broken: C-Cl",
    )
    bad = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|Ar-Cl|Ar-SH -> Ar-SR | bond_formed: C-O | bond_broken: C-Cl",
    )
    assert ok["reaction_type"] == "C_S_Coupling"
    assert bad["reaction_type"] == "Unknown"


def test_sonogashira_requires_c_c_bond_formation_in_taxonomy_constraints() -> None:
    ok = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|Ar-Cl|Ar-Alkynyl_terminal -> Ar-Alkynyl | bond_formed: C-C | bond_broken: C-Cl",
    )
    bad = validate_detection_with_crk_key(
        initial_detection="Unknown",
        initial_confidence=0.0,
        reaction_key="|Ar-Cl|Ar-Alkynyl_terminal -> Ar-Alkynyl | bond_formed: C-N | bond_broken: C-Cl",
    )
    assert ok["reaction_type"] == "Sonogashira"
    assert bad["reaction_type"] == "Unknown"
