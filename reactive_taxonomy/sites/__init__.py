"""Reactive-site detector registry."""

from . import electrophilic_centers, leaving_groups, pronucleophiles, transfer_groups

DETECTORS = {
    "leaving_group": leaving_groups.detect,
    "pronucleophile_XH": pronucleophiles.detect,
    "transfer_group": transfer_groups.detect,
    "electrophilic_center": electrophilic_centers.detect,
}

__all__ = ["DETECTORS"]
