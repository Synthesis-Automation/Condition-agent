"""
Rule-based electron flow templates keyed by mechanism archetype.

These heuristics are lightweight and deterministic, giving the agent useful
arrow-pushing descriptors without needing full quantum calculations.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional


_ELECTRON_FLOW_RULES: Dict[str, List[Dict[str, str]]] = {
    "sn2": [
        {
            "from": "nucleophile lone pair",
            "to": "sigma* (C-LG)",
            "description": "Backside attack forms the C-Nu bond as the C-LG bond breaks.",
        },
        {
            "from": "C-LG bond",
            "to": "leaving group",
            "description": "Leaving group departs with the electron pair.",
        },
    ],
    "sn1": [
        {
            "from": "C-LG bond",
            "to": "leaving group",
            "description": "Ionization generates a carbocation.",
        },
        {
            "from": "nucleophile lone pair",
            "to": "carbocation center",
            "description": "Nucleophile captures the planar carbocation.",
        },
    ],
    "addition_elimination_aromatic": [
        {
            "from": "nucleophile lone pair",
            "to": "ipso carbon",
            "description": "Nucleophile adds to the activated aromatic carbon (Meisenheimer intermediate).",
        },
        {
            "from": "Meisenheimer intermediate",
            "to": "leaving group",
            "description": "Aromaticity restored as LG departs.",
        },
    ],
    "oxidative_addition_reductive_elimination": [
        {
            "from": "metal d electrons",
            "to": "sigma* (C-X)",
            "description": "Pd/Cu inserts into the electrophile bond (oxidative addition).",
        },
        {
            "from": "nucleophile lone pair / transmetallated bond",
            "to": "metal center",
            "description": "Transmetalation or ligand exchange delivers Nu to the metal.",
        },
        {
            "from": "M-R / M-Nu bonds",
            "to": "product bond",
            "description": "Reductive elimination couples R and Nu, regenerating the catalyst.",
        },
    ],
    "transmetalation_coupling": [
        {
            "from": "metal d electrons",
            "to": "sigma* (C-X)",
            "description": "Oxidative addition of the electrophile.",
        },
        {
            "from": "M-C (organo partner)",
            "to": "Pd(II) center",
            "description": "Transmetalation transfers the nucleophilic fragment.",
        },
        {
            "from": "Pd-C bonds",
            "to": "biaryl bond",
            "description": "Reductive elimination forms the coupled product.",
        },
    ],
    "nucleophilic_acyl_substitution": [
        {
            "from": "nucleophile lone pair",
            "to": "carbonyl carbon",
            "description": "Addition into the carbonyl forms the tetrahedral intermediate.",
        },
        {
            "from": "tetrahedral intermediate",
            "to": "leaving group",
            "description": "Collapse expels the leaving group and restores the carbonyl.",
        },
    ],
    "beta_elimination": [
        {
            "from": "base lone pair",
            "to": "beta proton",
            "description": "Base abstracts beta-H anti to the leaving group.",
        },
        {
            "from": "C-H bond",
            "to": "C-C pi bond",
            "description": "Electrons flow to form the alkene as LG departs.",
        },
    ],
    "carbocation_elimination": [
        {
            "from": "C-LG bond",
            "to": "leaving group",
            "description": "Ionization yields the carbocation.",
        },
        {
            "from": "beta C-H bond",
            "to": "pi bond",
            "description": "Beta deprotonation forms the alkene.",
        },
    ],
    "pericyclic_cycloaddition": [
        {
            "from": "diene pi system",
            "to": "dienophile pi*",
            "description": "Concerted cyclic electron flow forms two sigma bonds simultaneously.",
        },
    ],
    "radical_chain": [
        {
            "from": "initiator bond",
            "to": "radical fragments",
            "description": "Initiation homolysis generates radicals.",
        },
        {
            "from": "radical center",
            "to": "substrate bond",
            "description": "Propagation forms product radical while regenerating the chain carrier.",
        },
    ],
    "photoredox_single_electron": [
        {
            "from": "photoexcited catalyst",
            "to": "substrate",
            "description": "SET step forms radical cation/anion.",
        },
        {
            "from": "radical intermediate",
            "to": "product bond",
            "description": "Radical capture/combination yields product.",
        },
    ],
    "copper_click_cycloaddition": [
        {
            "from": "Cu-acetylide",
            "to": "azide terminal nitrogen",
            "description": "Cu coordinates the alkyne and activates azide for cycloaddition.",
        },
        {
            "from": "dipolar intermediate",
            "to": "triazole ring",
            "description": "Ring closure delivers the 1,2,3-triazole and regenerates Cu(I).",
        },
    ],
}


def predict_electron_flow(
    mechanism_type: Optional[str],
    bond_changes: Optional[Dict[str, Any]] = None,
    electron_balance: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Return electron flow descriptors for the given mechanism type.
    """

    mechanism_key = (mechanism_type or "unknown").lower()
    arrows = list(_ELECTRON_FLOW_RULES.get(mechanism_key, []))

    notes: List[str] = []
    if not arrows:
        notes.append("No electron flow template for mechanism.")
    if bond_changes:
        broken = len(bond_changes.get("broken_bonds") or [])
        formed = len(bond_changes.get("formed_bonds") or [])
        notes.append(f"Bond evidence: {broken} broken / {formed} formed.")
    if electron_balance:
        atom_delta = electron_balance.get("atom_balance") or []
        if atom_delta:
            delta_str = ", ".join(
                f"{entry['atom']}: {entry['delta_lone_pairs']:+}"
                for entry in atom_delta[:5]
            )
            notes.append(f"Electron balance (Δ lone pairs): {delta_str}.")

    if arrows:
        diagram = " | ".join(
            f"{arrow.get('from', '?')} ==> {arrow.get('to', '?')}" for arrow in arrows
        )
    else:
        diagram = ""

    return {
        "mechanism_type": mechanism_key,
        "arrows": arrows,
        "notes": notes,
        "diagram": diagram,
    }


__all__ = ["predict_electron_flow"]
