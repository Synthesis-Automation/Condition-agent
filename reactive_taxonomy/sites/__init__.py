"""Reactive-site detector registry."""

from . import aromatic_ch, dipolar_groups, electrophilic_centers, leaving_groups, pronucleophiles, transfer_groups, unsaturated_bonds

DETECTORS = {
    "leaving_group": leaving_groups.detect,
    "pronucleophile_XH": pronucleophiles.detect,
    "transfer_group": transfer_groups.detect,
    "electrophilic_center": electrophilic_centers.detect,
    "aromatic_CH": aromatic_ch.detect,
    "unsaturated_bond": unsaturated_bonds.detect,
    "dipolar_group": dipolar_groups.detect,
}

__all__ = ["DETECTORS"]
