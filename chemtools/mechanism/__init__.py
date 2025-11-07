"""
Lightweight, rule-based mechanism heuristics for the agent mechanism tool.

Option-1 scope: map detected reaction families onto canonical mechanism
archetypes, surface informative descriptions, and produce a structured
step list that downstream agents can narrate. This module intentionally
avoids heavyweight computation so it is always available.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional


_DETAIL_LEVELS = {"summary", "standard", "high"}

_FAMILY_SYNONYMS: Dict[str, str] = {
    "buchwald_cn": "buchwald_hartwig_c_n",
    "buchwald-cn": "buchwald_hartwig_c_n",
    "buchwaldhartwig": "buchwald_hartwig_c_n",
    "buchwald_hartwig": "buchwald_hartwig_c_n",
    "buchwaldhartwig_cn": "buchwald_hartwig_c_n",
    "c_n_coupling": "cn_coupling",
    "c-n_coupling": "cn_coupling",
    "cu_cn_coupling": "cn_coupling",
    "cncoupling": "cn_coupling",
    "sn2_alkyl": "sn2_alkylation",
    "sn2alkylation": "sn2_alkylation",
    "williamsonether": "williamson_ether",
    "amide_coupling": "amide_formation",
    "amidation": "amide_formation",
    "sn_ar": "snar",
    "s_nar": "snar",
    "aromatic_nucleophilic_substitution": "snar",
    "nucleophilic_aromatic_substitution": "snar",
    "pericyclic": "diels_alder",
    "diels-alder": "diels_alder",
    "diels_alder_reaction": "diels_alder",
    "radical_addition": "radical_chain",
    "radical_polymerization": "radical_chain",
    "photoredox_coupling": "photoredox_single_electron",
    "click_reaction": "cu_azide_alkyne_click",
}

_FAMILY_TO_MECHANISM: Dict[str, str] = {
    "cn_coupling": "oxidative_addition_reductive_elimination",
    "buchwald_hartwig_c_n": "oxidative_addition_reductive_elimination",
    "buchwald_hartwig": "oxidative_addition_reductive_elimination",
    "buchwald_hartwig_cn": "oxidative_addition_reductive_elimination",
    "ullmann": "oxidative_addition_reductive_elimination",
    "suzuki": "transmetalation_coupling",
    "suzuki_miyaura": "transmetalation_coupling",
    "suzuki_coupling": "transmetalation_coupling",
    "negishi": "transmetalation_coupling",
    "stille": "transmetalation_coupling",
    "sonogashira": "transmetalation_coupling",
    "amide_formation": "nucleophilic_acyl_substitution",
    "nucleophilic_acyl_substitution": "nucleophilic_acyl_substitution",
    "knoevenagel": "nucleophilic_acyl_substitution",
    "sn2_alkylation": "sn2",
    "williamson_ether": "sn2",
    "sn2": "sn2",
    "sn1": "sn1",
    "snar": "addition_elimination_aromatic",
    "e1": "carbocation_elimination",
    "e2": "beta_elimination",
    "elimination": "beta_elimination",
    "heck": "migratory_insertion_beta_hydride",
    "diels_alder": "pericyclic_cycloaddition",
    "cope": "pericyclic_cycloaddition",
    "radical_chain": "radical_chain",
    "photoredox_single_electron": "photoredox_single_electron",
    "cu_azide_alkyne_click": "copper_click_cycloaddition",
}

_MECHANISM_DESCRIPTIONS: Dict[str, str] = {
    "oxidative_addition_reductive_elimination": (
        "Pd/Cu-catalyzed cross-coupling: oxidative addition of the electrophile, "
        "transmetalation or ligand exchange, then reductive elimination to form the bond."
    ),
    "transmetalation_coupling": (
        "Cross-coupling via oxidative addition, transmetalation from an organometallic "
        "partner, and reductive elimination to deliver the coupled product."
    ),
    "nucleophilic_acyl_substitution": (
        "Carbonyl addition-elimination: nucleophile attacks the carbonyl, a "
        "tetrahedral intermediate forms, then the leaving group departs."
    ),
    "sn2": (
        "Concerted backside attack where the nucleophile displaces the leaving group "
        "with inversion of configuration."
    ),
    "sn1": (
        "Two-step substitution via carbocation formation followed by nucleophile capture."
    ),
    "beta_elimination": (
        "Base-induced removal of a proton beta to a leaving group, forming a new pi bond."
    ),
    "carbocation_elimination": (
        "E1 elimination: leaving group departure to form a carbocation, then beta deprotonation yields the alkene."
    ),
    "addition_elimination_aromatic": (
        "SNAr addition-elimination: nucleophile adds to the activated aromatic ring, "
        "forming a Meisenheimer intermediate before expelling the leaving group."
    ),
    "migratory_insertion_beta_hydride": (
        "Heck-type sequence: oxidative addition, alkene migratory insertion, "
        "beta-hydride elimination, and catalyst regeneration."
    ),
    "pericyclic_cycloaddition": (
        "Concerted pericyclic pathway (e.g., Diels-Alder) where bonds form and break "
        "in a single transition state guided by orbital symmetry."
    ),
    "radical_chain": (
        "Radical chain: initiation to generate radicals, propagation to build product, "
        "and termination via radical recombination or quenching."
    ),
    "photoredox_single_electron": (
        "Photoredox single-electron transfer: photocatalyst excitation, SET events, "
        "radical or cationic intermediates, and catalyst turnover."
    ),
    "copper_click_cycloaddition": (
        "CuAAC click reaction: copper acetylide formation, azide activation, "
        "and 1,3-dipolar cycloaddition to a triazole."
    ),
}

_MECHANISM_STEP_LIBRARY: Dict[str, List[Dict[str, Any]]] = {
    "oxidative_addition_reductive_elimination": [
        {
            "title": "Oxidative addition",
            "description": "Pd(0)/Cu(I) inserts into the carbon-leaving-group bond.",
            "electron_flow": "LnM + R-X -> LnM(II)(R)(X)",
        },
        {
            "title": "Transmetalation / ligand exchange",
            "description": "Organometallic partner or amide nucleophile replaces X on the metal.",
            "electron_flow": "LnM(II)(R)(X) + Nu -> LnM(II)(R)(Nu) + X-",
        },
        {
            "title": "Reductive elimination",
            "description": "Coupling partners form the new bond; catalyst returns to lower oxidation state.",
            "electron_flow": "LnM(II)(R)(Nu) -> LnM(0) + R-Nu",
        },
    ],
    "transmetalation_coupling": [
        {
            "title": "Oxidative addition",
            "description": "Pd(0) inserts into the C-X bond of the electrophile.",
            "electron_flow": "Pd(0) + Ar-X -> Pd(II)(Ar)(X)",
        },
        {
            "title": "Transmetalation",
            "description": "Organoboron/zinc transfers the nucleophilic fragment to Pd.",
            "electron_flow": "Pd(II)(Ar)(X) + Ar'-B(OH)2 -> Pd(II)(Ar)(Ar') + X-B(OH)2",
        },
        {
            "title": "Reductive elimination",
            "description": "Biaryl bond forms and Pd(0) is regenerated.",
            "electron_flow": "Pd(II)(Ar)(Ar') -> Pd(0) + Ar-Ar'",
        },
    ],
    "nucleophilic_acyl_substitution": [
        {
            "title": "Nucleophilic attack",
            "description": "Nucleophile adds to the electrophilic carbonyl carbon.",
            "electron_flow": "Nu: -> C=O, forming a tetrahedral intermediate.",
        },
        {
            "title": "Tetrahedral intermediate collapse",
            "description": "Leaving group departs, restoring the carbonyl.",
            "electron_flow": "Intermediate -> Product + LG-",
        },
    ],
    "sn2": [
        {
            "title": "Backside attack",
            "description": "Nucleophile approaches antiperiplanar to the leaving group.",
            "electron_flow": "Nu: -> C-LG, bond forms as C-LG breaks.",
        },
    ],
    "sn1": [
        {
            "title": "Carbocation formation",
            "description": "C-LG bond ionizes to yield a carbocation.",
            "electron_flow": "C-LG -> C+ + LG-",
        },
        {
            "title": "Nucleophile capture",
            "description": "Nucleophile attacks the planar carbocation.",
            "electron_flow": "Nu: -> C+, forming C-Nu.",
        },
    ],
    "carbocation_elimination": [
        {
            "title": "Leaving group departure",
            "description": "Beta-heteroatom bond cleaves, generating a carbocation.",
        },
        {
            "title": "Beta deprotonation",
            "description": "Base removes a beta proton to furnish the alkene.",
        },
    ],
    "beta_elimination": [
        {
            "title": "Beta-proton abstraction",
            "description": "Base removes a beta-hydrogen to form a C=C bond.",
            "electron_flow": "Base: -> Hbeta; electrons flow to form pi bond as LG departs.",
        },
    ],
    "addition_elimination_aromatic": [
        {
            "title": "Nucleophilic addition",
            "description": "Nucleophile adds to the activated aromatic carbon bearing the leaving group.",
        },
        {
            "title": "Meisenheimer intermediate collapse",
            "description": "Aromaticity is restored as the leaving group departs.",
        },
    ],
    "migratory_insertion_beta_hydride": [
        {
            "title": "Oxidative addition",
            "description": "Pd(0) adds into the C-X bond of the vinyl/aryl halide.",
        },
        {
            "title": "Migratory insertion",
            "description": "Coordinated alkene inserts into the Pd-C bond.",
        },
        {
            "title": "Beta-hydride elimination",
            "description": "Elimination delivers the substituted alkene and regenerates Pd(0).",
        },
    ],
    "pericyclic_cycloaddition": [
        {
            "title": "Concerted cycloaddition",
            "description": "Diene and dienophile orbitals reorganize in a single step to make two sigma bonds.",
        },
    ],
    "radical_chain": [
        {
            "title": "Initiation",
            "description": "Heat, light, or peroxide forms the first radical.",
        },
        {
            "title": "Propagation",
            "description": "Radicals add/abstract to build product and regenerate the radical.",
        },
        {
            "title": "Termination",
            "description": "Radicals recombine or are quenched, ending the chain.",
        },
    ],
    "photoredox_single_electron": [
        {
            "title": "Photocatalyst excitation",
            "description": "Catalyst absorbs light to reach an excited state.",
        },
        {
            "title": "Single-electron transfer",
            "description": "Electron transfer forms radical or cationic intermediates.",
        },
        {
            "title": "Product-forming steps",
            "description": "Radical capture, recombination, or elimination delivers the product.",
        },
    ],
    "copper_click_cycloaddition": [
        {
            "title": "Copper acetylide formation",
            "description": "Cu(I) coordinates and deprotonates the terminal alkyne.",
        },
        {
            "title": "Azide activation and cycloaddition",
            "description": "Copper-bound azide undergoes 1,3-dipolar cycloaddition to the acetylide.",
        },
    ],
}


def classify_mechanism_simple(
    reaction_family: Optional[str],
    bond_changes: Optional[Dict[str, Any]] = None,
    detail_level: str = "standard",
) -> Dict[str, Any]:
    """
    Classify a mechanism archetype from a detected reaction family and bond data.

    Args:
        reaction_family: Canonical family label from chemtools.detect_reaction.
        bond_changes: Optional bond analysis payload for evidence.
        detail_level: "summary", "standard", or "high".

    Returns:
        Dict with mechanism_type, confidence, description, steps, warnings,
        and evidence references suitable for agent consumption.
    """

    normalized_level = detail_level if detail_level in _DETAIL_LEVELS else "standard"
    normalized_family = _normalize_family(reaction_family)
    mechanism_type = _FAMILY_TO_MECHANISM.get(
        normalized_family,
        "unknown",
    )

    confidence = 0.75 if mechanism_type != "unknown" else 0.35
    description = _MECHANISM_DESCRIPTIONS.get(
        mechanism_type,
        "Mechanism details not available for this reaction family.",
    )

    steps = _render_steps(mechanism_type, normalized_level)
    bond_step = _summarize_bond_changes(bond_changes)
    if bond_step:
        steps.append(bond_step)

    warnings: List[str] = []
    if mechanism_type == "unknown":
        warnings.append("Mechanism type could not be inferred from the detected family.")
    if bond_changes and not bond_step:
        warnings.append("Bond analysis returned no clear bond changes.")

    evidence_refs: List[Dict[str, Any]] = []
    if reaction_family:
        evidence_refs.append(
            {
                "source": "reaction_detection",
                "detail": {"family": reaction_family},
            }
        )
    if bond_changes:
        evidence_refs.append(
            {
                "source": "bond_analysis",
                "detail": {
                    "broken_bonds": bond_changes.get("broken_bonds", []),
                    "formed_bonds": bond_changes.get("formed_bonds", []),
                    "changed_atoms": list(bond_changes.get("changed_atoms", []))
                    if isinstance(bond_changes.get("changed_atoms"), (list, set, tuple))
                    else bond_changes.get("changed_atoms"),
                },
            }
        )

    return {
        "mechanism_type": mechanism_type,
        "description": description,
        "confidence": confidence,
        "detail_level": normalized_level,
        "steps": steps,
        "warnings": warnings,
        "evidence_refs": evidence_refs,
    }


def _normalize_family(family: Optional[str]) -> str:
    if not family:
        return "unknown"
    normalized = family.strip().lower().replace("-", "_").replace(" ", "_")
    normalized = normalized.replace("__", "_")
    if normalized in _FAMILY_SYNONYMS:
        normalized = _FAMILY_SYNONYMS[normalized]
    return normalized


def _render_steps(mechanism_type: str, detail_level: str) -> List[Dict[str, Any]]:
    library = _MECHANISM_STEP_LIBRARY.get(mechanism_type, [])
    if not library:
        return []
    limit = 1 if detail_level == "summary" else len(library)
    if detail_level == "high":
        # Duplicate list to avoid mutating constants and annotate level.
        steps = [
            {**step, "detail_level": "high"}
            for step in library[:limit]
        ]
    else:
        steps = [dict(step) for step in library[:limit]]
    return steps


def _summarize_bond_changes(
    bond_changes: Optional[Dict[str, Any]]
) -> Optional[Dict[str, Any]]:
    if not bond_changes:
        return None
    broken = bond_changes.get("broken_bonds") or []
    formed = bond_changes.get("formed_bonds") or []
    changed_atoms = bond_changes.get("changed_atoms") or []
    if not broken and not formed and not changed_atoms:
        return None
    description = (
        f"{len(broken)} bond(s) break; {len(formed)} bond(s) form."
    )
    return {
        "title": "Bond change summary",
        "description": description,
        "electron_flow": "Derived from RXNMapper/MCS analysis.",
        "key_atoms": list(changed_atoms) if isinstance(changed_atoms, (list, tuple, set)) else changed_atoms,
    }


from .electron_flow import predict_electron_flow
from .intermediates import predict_intermediates

__all__ = ["classify_mechanism_simple", "predict_electron_flow", "predict_intermediates"]
