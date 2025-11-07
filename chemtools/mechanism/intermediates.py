"""
Heuristic intermediate prediction aligned with mechanism archetypes.
"""

from __future__ import annotations

from typing import Dict, List, Optional


_INTERMEDIATE_RULES: Dict[str, List[Dict[str, str]]] = {
    "sn1": [
        {
            "name": "carbocation",
            "description": "Planar carbocation formed after leaving group departure.",
            "stability": "Stabilized by resonance/hyperconjugation.",
        }
    ],
    "sn2": [
        {
            "name": "pentavalent transition structure",
            "description": "Single concerted step; no discrete intermediate, but backside TS noted.",
            "stability": "Lower barrier with good leaving group & polar aprotic solvent.",
        }
    ],
    "oxidative_addition_reductive_elimination": [
        {
            "name": "Pd(II) oxidative addition complex",
            "description": "LnPd(II)(R)(X) formed after oxidative addition.",
            "stability": "Favored by electron-rich ligands.",
        },
        {
            "name": "transmetalated Pd(II)",
            "description": "LnPd(II)(R)(Nu) prior to reductive elimination.",
            "stability": "Stability depends on ligand sterics/electronics.",
        },
    ],
    "transmetalation_coupling": [
        {
            "name": "Pd(II) diorganometallic",
            "description": "Both coupling partners bound before reductive elimination.",
            "stability": "Stabilized by phosphine or NHC ligands.",
        }
    ],
    "nucleophilic_acyl_substitution": [
        {
            "name": "tetrahedral intermediate",
            "description": "Addition to carbonyl generates sp3 carbon bearing both nucleophile and LG.",
            "stability": "Enhanced by electron-withdrawing acyl groups.",
        }
    ],
    "addition_elimination_aromatic": [
        {
            "name": "Meisenheimer complex",
            "description": "Anionic sigma complex with disrupted aromaticity.",
            "stability": "Needs strong EWG ortho/para to LG.",
        }
    ],
    "beta_elimination": [
        {
            "name": "anti-periplanar arrangement",
            "description": "Pre-elimination conformer with beta-H anti to leaving group.",
            "stability": "Conformational preference controls rate.",
        }
    ],
    "pericyclic_cycloaddition": [
        {
            "name": "concerted transition state",
            "description": "Cyclic TS with synchronous bond formation.",
            "stability": "Thermal vs photochemical control via Woodward-Hoffmann rules.",
        }
    ],
    "radical_chain": [
        {
            "name": "radical propagation intermediates",
            "description": "Substrate-centered radicals after hydrogen abstraction or addition.",
            "stability": "Benzylic/allylic radicals favored.",
        }
    ],
    "photoredox_single_electron": [
        {
            "name": "radical cation/anion",
            "description": "Generated via single-electron transfer with photocatalyst.",
            "stability": "Depends on redox potential and solvent.",
        }
    ],
    "copper_click_cycloaddition": [
        {
            "name": "copper-acetylide",
            "description": "Cu(I) bound to deprotonated terminal alkyne.",
            "stability": "Supported by ligands (TBTA, bipyridine).",
        },
        {
            "name": "metallacyclic azide complex",
            "description": "Cu-coordinated 1,3-dipole prior to ring closure.",
            "stability": "Short-lived, collapses to triazole.",
        },
    ],
}


def predict_intermediates(
    mechanism_type: Optional[str],
    reaction_family: Optional[str] = None,
    context: Optional[Dict[str, str]] = None,
) -> List[Dict[str, str]]:
    """
    Return representative intermediates for the mechanism type.
    """

    mechanism_key = (mechanism_type or "unknown").lower()
    intermediates = list(_INTERMEDIATE_RULES.get(mechanism_key, []))

    if not intermediates and reaction_family:
        # fallback: mention unknown intermediate yet note family
        intermediates = [
            {
                "name": f"{reaction_family} intermediate",
                "description": "No curated intermediate available; review reaction context.",
                "stability": "Unknown",
            }
        ]

    context = context or {}
    if intermediates and context.get("catalyst"):
        for entry in intermediates:
            entry.setdefault("notes", f"Consistent with catalyst {context['catalyst']}.")

    return intermediates


__all__ = ["predict_intermediates"]
