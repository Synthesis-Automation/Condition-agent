from __future__ import annotations

from functools import lru_cache
from typing import Iterable, Sequence, Any, List

import numpy as np

from chemtools.integrations import molpipeline as mp

__all__ = ["morgan_fingerprint", "physchem_features"]


def _ensure_smiles_sequence(smiles_like: Iterable[Any] | Sequence[str] | str) -> List[str]:
    """Normalize caller input to a list of SMILES strings."""

    if isinstance(smiles_like, str):
        seq = [smiles_like]
    elif isinstance(smiles_like, Sequence):
        seq = list(smiles_like)
    else:
        seq = list(smiles_like)
    if not seq:
        raise ValueError("No SMILES strings provided.")
    return seq


@lru_cache(maxsize=16)
def _cached_morgan_pipeline(n_bits: int, radius: int, n_jobs: int):
    return mp.build_morgan_pipeline(n_bits=n_bits, radius=radius, n_jobs=n_jobs)


@lru_cache(maxsize=16)
def _cached_physchem_pipeline(descriptor_key: tuple[str, ...] | None, n_jobs: int):
    descriptor_list = list(descriptor_key) if descriptor_key is not None else None
    return mp.build_physchem_pipeline(descriptor_list, n_jobs=n_jobs)


def morgan_fingerprint(
    smiles_like: Iterable[Any] | Sequence[str] | str,
    *,
    n_bits: int = 2048,
    radius: int = 2,
    n_jobs: int = 1,
    return_sparse: bool = True,
):
    """Generate Morgan fingerprints for the provided molecules via MolPipeline."""

    if not mp.is_available():
        raise RuntimeError("MolPipeline is not available; install the molpipeline extra first.")
    smiles_seq = _ensure_smiles_sequence(smiles_like)
    pipeline = _cached_morgan_pipeline(n_bits, radius, n_jobs)
    matrix = mp.transform_smiles(smiles_seq, pipeline)
    if return_sparse:
        return matrix
    if hasattr(matrix, "toarray"):
        dense = matrix.toarray()
    else:
        dense = matrix
    return np.asarray(dense, dtype=float, copy=True)


def physchem_features(
    smiles_like: Iterable[Any] | Sequence[str] | str,
    *,
    descriptor_list: Sequence[str] | None = None,
    n_jobs: int = 1,
) -> np.ndarray:
    """Return RDKit physicochemical descriptors using MolPipeline convenience wrappers."""

    if not mp.is_available():
        raise RuntimeError("MolPipeline is not available; install the molpipeline extra first.")
    smiles_seq = _ensure_smiles_sequence(smiles_like)
    descriptor_key = tuple(descriptor_list) if descriptor_list is not None else None
    pipeline = _cached_physchem_pipeline(descriptor_key, n_jobs)
    matrix = mp.transform_smiles(smiles_seq, pipeline)
    if hasattr(matrix, "toarray"):
        return np.asarray(matrix.toarray(), dtype=float)
    return np.asarray(matrix, dtype=float)
