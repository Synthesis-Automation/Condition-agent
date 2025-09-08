from __future__ import annotations

"""
Lightweight DRFP-based reaction similarity utilities.

This module is optional: if the `drfp` package (and NumPy) are not installed,
all functions degrade gracefully and return None or 0.0 so callers can keep
working with categorical-only similarity.

Public functions:
- drfp_available() -> bool
- encode_drfp(rsmi: str, n_bits: int = 4096, radius: int = 3) -> Any | None
- tanimoto(a, b) -> float
"""

from typing import Any, Optional
from functools import lru_cache


def _import_drfp():
    try:
        from drfp import DrfpEncoder  # type: ignore
        return DrfpEncoder
    except Exception:
        return None


def _import_numpy():
    try:
        import numpy as np  # type: ignore
        return np
    except Exception:
        return None


_PRECOMP: dict | None = {}


def drfp_available() -> bool:
    return _import_drfp() is not None and _import_numpy() is not None


def encode_drfp(rsmi: str, n_bits: int = 4096, radius: int = 3) -> Optional[Any]:
    """Encode a reaction SMILES with DRFP.

    Returns a NumPy array of dtype uint8 with 0/1 bits, or None if unavailable.
    """
    DrfpEncoder = _import_drfp()
    np = _import_numpy()
    if DrfpEncoder is None or np is None:
        return None
    try:
        res = DrfpEncoder.encode(
            str(rsmi),
            n_folded_length=int(n_bits),
            min_radius=0,
            radius=int(radius),
            rings=True,
            mapping=False,
            atom_index_mapping=False,
            root_central_atom=True,
            include_hydrogens=False,
            show_progress_bar=False,
        )
        # drfp>=0.3 returns a list for iterable input; for string input it may return a single array
        try:
            import numpy as _np  # type: ignore
        except Exception:
            _np = None
        if isinstance(res, list):
            arr = res[0] if res else None
        else:
            arr = res
        if arr is None:
            return None
        if _np is not None:
            arr = _np.asarray(arr).astype('uint8')
        return arr
    except Exception:
        return None


def tanimoto(a: Any, b: Any) -> float:
    """Compute Tanimoto similarity for two DRFP bit arrays.

    Expects NumPy arrays with 0/1 values. Returns 0.0 if inputs are invalid.
    """
    np = _import_numpy()
    if np is None:
        return 0.0
    try:
        if a is None or b is None:
            return 0.0
        # Ensure boolean-ish arrays
        a = (a > 0).astype('uint8')
        b = (b > 0).astype('uint8')
        inter = int(np.sum((a & b) == 1))
        denom = int(np.sum(a == 1) + np.sum(b == 1) - inter)
        return 0.0 if denom <= 0 else float(inter) / float(denom)
    except Exception:
        return 0.0


@lru_cache(maxsize=200000)
def encode_drfp_cached(rsmi: str, n_bits: int = 4096, radius: int = 3) -> Optional[Any]:
    """LRU-cached DRFP encoding by reaction SMILES and parameters.

    Returns a NumPy array or None. Cache size is generous to accommodate
    medium datasets. Arrays are kept in-memory for fast reuse.
    """
    # Prefer precomputed vectors when available
    try:
        key = (str(rsmi), int(n_bits), int(radius))
        if isinstance(_PRECOMP, dict) and key in _PRECOMP:
            return _PRECOMP[key]  # type: ignore[index]
    except Exception:
        pass
    return encode_drfp(rsmi, n_bits=n_bits, radius=radius)


def load_precomputed_npz(path: str) -> dict:
    """Load a precomputed NPZ bundle into the in-memory store.

    Expects arrays: 'fps' (N x n_bits, uint8), 'keys' (N strings). Optional scalars
    'n_bits', 'radius'. Returns a short status dict.
    """
    np = _import_numpy()
    if np is None:
        return {"ok": False, "error": "numpy not available"}
    try:
        data = np.load(path, allow_pickle=True)
    except Exception as e:
        return {"ok": False, "error": f"load failed: {e}"}
    try:
        fps = data["fps"]
        keys = data["keys"]
    except Exception as e:
        return {"ok": False, "error": f"missing arrays: {e}"}
    try:
        n_bits = int(data["n_bits"].item()) if "n_bits" in data else int(fps.shape[1])
    except Exception:
        n_bits = int(fps.shape[1])
    try:
        radius = int(data["radius"].item()) if "radius" in data else 3
    except Exception:
        radius = 3
    loaded = 0
    global _PRECOMP
    if _PRECOMP is None:
        _PRECOMP = {}
    for i in range(len(keys)):
        try:
            rsmi = str(keys[i])
            vec = fps[i]
            _PRECOMP[(rsmi, n_bits, radius)] = vec  # type: ignore[index]
            loaded += 1
        except Exception:
            continue
    return {"ok": True, "count": loaded, "n_bits": n_bits, "radius": radius}
