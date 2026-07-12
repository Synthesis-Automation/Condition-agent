"""Small, dependency-local RDKit helpers.

The standalone taxonomy depends directly on RDKit and never imports the legacy
``chemtools`` package. Parsing failures are intentionally returned as ``None``
because invalid user input is represented by the public analysis models.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator, Optional


def _chem_module() -> Any:
    try:
        from rdkit import Chem  # type: ignore
    except ImportError:
        return None
    return Chem


@contextmanager
def _suppress_parse_logging() -> Iterator[None]:
    try:
        from rdkit import rdBase  # type: ignore
    except ImportError:
        yield
        return
    blocker = getattr(rdBase, "BlockLogs", None)
    if blocker is None:
        yield
        return
    with blocker():
        yield


def rdkit_available() -> bool:
    """Return whether RDKit can be imported."""
    return _chem_module() is not None


def parse_smiles(smiles: str) -> Any:
    """Parse a SMILES string, returning ``None`` for invalid input."""
    chem = _chem_module()
    if chem is None:
        return None
    try:
        with _suppress_parse_logging():
            return chem.MolFromSmiles(str(smiles or ""))
    except Exception:
        return None


def mol_to_canonical_smiles(mol: Any) -> Optional[str]:
    """Serialize an RDKit molecule to canonical isomeric SMILES."""
    chem = _chem_module()
    if chem is None or mol is None:
        return None
    try:
        return str(chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    except Exception:
        return None
