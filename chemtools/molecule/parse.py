"""Molecule parsing helpers."""

from __future__ import annotations

from chemtools.core.rdkit import canonical_smiles, parse_smiles, rdkit_available

from .models import MoleculeParseResult


def parse_molecule(smiles: str) -> MoleculeParseResult:
    """Parse a molecule SMILES string into a small domain result."""
    text = str(smiles or "").strip()
    available = rdkit_available()
    if not text:
        return MoleculeParseResult(
            input_smiles=text,
            canonical_smiles=None,
            rdkit_available=available,
            valid=False,
            error="EMPTY_SMILES",
        )
    if not available:
        return MoleculeParseResult(
            input_smiles=text,
            canonical_smiles=text,
            rdkit_available=False,
            valid=True,
        )
    mol = parse_smiles(text)
    if mol is None:
        return MoleculeParseResult(
            input_smiles=text,
            canonical_smiles=None,
            rdkit_available=True,
            valid=False,
            error="INVALID_SMILES",
        )
    return MoleculeParseResult(
        input_smiles=text,
        canonical_smiles=canonical_smiles(text),
        rdkit_available=True,
        valid=True,
    )


__all__ = ["parse_molecule"]
