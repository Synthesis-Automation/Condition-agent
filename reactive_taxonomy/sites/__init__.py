"""Reactive-site detector registry."""

from . import (
    addition_donors,
    anionic_nucleophiles,
    aromatic_ch,
    dipolar_groups,
    electrophilic_centers,
    eliminable_pairs,
    heteroatom_bonds,
    leaving_groups,
    pronucleophiles,
    transfer_groups,
    unsaturated_bonds,
)

DETECTORS = {
    "leaving_group": leaving_groups.detect,
    "pronucleophile_XH": pronucleophiles.detect,
    "nucleophile_anion": anionic_nucleophiles.detect,
    "transfer_group": transfer_groups.detect,
    "addition_donor": addition_donors.detect,
    "eliminable_pair": eliminable_pairs.detect,
    "electrophilic_center": electrophilic_centers.detect,
    "aromatic_CH": aromatic_ch.detect,
    "unsaturated_bond": unsaturated_bonds.detect,
    "dipolar_group": dipolar_groups.detect,
    "heteroatom_bond": heteroatom_bonds.detect,
}

__all__ = ["DETECTORS"]
