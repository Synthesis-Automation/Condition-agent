"""Small, dependency-local RDKit helpers.

The standalone taxonomy depends directly on RDKit. Parsing failures are
intentionally returned as ``None`` because invalid user input is represented by
the public analysis models.
"""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from typing import Any, Iterable, Iterator, Optional


@dataclass
class MoleculeParseCache:
    """Invocation-owned immutable binary parses; every consumer gets a fresh Mol."""

    max_entries: int = 512
    entries: dict[str, bytes | None] = field(default_factory=dict)


_PARSE_SCOPE: ContextVar[MoleculeParseCache | None] = ContextVar("taxonomy_parse_scope", default=None)


@contextmanager
def molecule_parse_scope(cache: MoleculeParseCache) -> Iterator[None]:
    """Reuse parsing within an explicit invocation, isolated by context."""
    token = _PARSE_SCOPE.set(cache)
    try:
        yield
    finally:
        _PARSE_SCOPE.reset(token)


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
        key = str(smiles or "")
        cache = _PARSE_SCOPE.get()
        if cache is not None and key in cache.entries:
            stored = cache.entries[key]
            return chem.Mol(stored) if stored is not None else None
        with _suppress_parse_logging():
            molecule = chem.MolFromSmiles(key)
        if cache is not None:
            if len(cache.entries) >= cache.max_entries:
                cache.entries.pop(next(iter(cache.entries)))
            cache.entries[key] = molecule.ToBinary(chem.PropertyPickleOptions.AllProps) if molecule is not None else None
        return molecule
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


def prepare_fragment_serialization_copy(
    mol: Any,
    atom_indices: Iterable[int],
) -> Any:
    """Copy a molecule and clear stereo that crosses a fragment boundary.

    RDKit's fragment canonicalizer inspects double-bond stereo annotations on
    the parent molecule. A stereobond wholly outside ``atomsToUse`` can violate
    its traversal precondition. Stereo wholly inside the selected induced
    subgraph is retained; only annotations that the fragment cannot represent
    are cleared.
    """
    chem = _chem_module()
    if chem is None or mol is None:
        return None
    selected = {int(value) for value in atom_indices}
    copied = chem.Mol(mol)
    for bond in copied.GetBonds():
        if (
            int(bond.GetBeginAtomIdx()) in selected
            and int(bond.GetEndAtomIdx()) in selected
        ):
            continue
        bond.SetStereo(chem.BondStereo.STEREONONE)
        bond.SetBondDir(chem.BondDir.NONE)
    return copied
