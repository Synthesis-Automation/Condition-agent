"""
Molecular featurizer (general representation).

Implements a stable entry point `featurize(electrophile, nucleophile)`
without emitting deprecation warnings (unlike the legacy `ullmann`
module). Internally reuses the cached core implementation from
`chemtools.featurizers.ullmann` and conditionally attaches role-aware
vectors when the optional `chem_feats` package is available.
"""

from __future__ import annotations

from typing import Dict, Any
import os

from .ullmann import _featurize_cached as _u_cached  # core implementation

try:  # Optional role-aware vectors
    from chem_feats import featurize_mol as _role_feat  # type: ignore
    _HAS_ROLE_FEATS = True
except Exception:  # pragma: no cover
    _role_feat = None  # type: ignore
    _HAS_ROLE_FEATS = False


def featurize(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    """Featurize electrophile/nucleophile pair without deprecation noise.

    Returns a dict including coarse LG/nucleophile classes and a `bin` key.
    When role-aware features are available, attaches a `role_aware` block with
    molecule-level vectors for the electrophile and nucleophile.
    """
    base = _u_cached(electrophile, nucleophile)
    out = dict(base)
    # Attach role-aware only when explicitly enabled via env (default off for speed)
    attach_flag = (os.environ.get("CHEMTOOLS_ATTACH_ROLE_AWARE", "").strip().lower() in {"1", "true", "yes", "on"})
    if attach_flag and _HAS_ROLE_FEATS and _role_feat is not None:
        try:
            elec_ra = _role_feat(electrophile, roles=["aryl_halide"])  # type: ignore[arg-type]
            nuc_ra = _role_feat(nucleophile, roles=["amine", "alcohol"])  # type: ignore[arg-type]
            out["role_aware"] = {
                "electrophile": elec_ra,
                "nucleophile": nuc_ra,
            }
        except Exception:
            pass
    return out


__all__ = ["featurize"]
