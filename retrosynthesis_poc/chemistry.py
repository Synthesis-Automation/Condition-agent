"""Small deterministic chemistry helpers for the retrosynthesis POC."""

from __future__ import annotations

import hashlib
from functools import lru_cache
from typing import Iterable, Optional, Tuple

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def digest(namespace: str, *values: str) -> str:
    """Return a deterministic namespaced SHA-256 identifier."""

    payload = "\0".join(values).encode("utf-8")
    return f"{namespace}:{hashlib.sha256(payload).hexdigest()}"


def molecule_without_maps(smiles: str) -> Optional[object]:
    """Parse SMILES and return a sanitized copy without atom-map numbers."""

    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        return None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        Chem.SanitizeMol(molecule)
    except Exception:
        return None
    return molecule


def canonical_smiles(smiles: str) -> Optional[str]:
    """Canonicalize one molecule or dot-separated molecule collection."""

    molecule = molecule_without_maps(smiles)
    if molecule is None:
        return None
    fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
    values = sorted(
        Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=True)
        for fragment in fragments
    )
    return ".".join(values)


def split_reaction_smiles(reaction_smiles: str) -> Optional[Tuple[str, str]]:
    """Return mapped reactant and product sides for a two-part reaction."""

    if reaction_smiles.count(">>") != 1:
        return None
    reactants, products = reaction_smiles.split(">>")
    if not reactants or not products:
        return None
    return reactants, products


def atom_map_numbers(smiles: str) -> set[int]:
    """Return positive atom-map numbers present in a SMILES collection."""

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return set()
    return {
        int(atom.GetAtomMapNum())
        for atom in molecule.GetAtoms()
        if int(atom.GetAtomMapNum()) > 0
    }


def contributing_reactants(
    reactant_smiles: str,
    product_smiles: str,
) -> Optional[str]:
    """Drop disconnected spectators that contribute no mapped product atom."""

    product_maps = atom_map_numbers(product_smiles)
    if not product_maps:
        return None
    values = []
    for component in reactant_smiles.split("."):
        if atom_map_numbers(component).intersection(product_maps):
            canonical = canonical_smiles(component)
            if canonical is None:
                return None
            values.append(canonical)
    if not values:
        return None
    return ".".join(sorted(values))


@lru_cache(maxsize=1)
def _morgan_generator() -> object:
    return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@lru_cache(maxsize=50_000)
def fingerprint(smiles: str) -> Optional[object]:
    """Return a cached Morgan fingerprint for canonical SMILES."""

    molecule = molecule_without_maps(smiles)
    if molecule is None:
        return None
    return _morgan_generator().GetFingerprint(molecule)


def maximum_similarity(query_smiles: str, references: Iterable[str]) -> float:
    """Return maximum Tanimoto similarity to any valid reference."""

    query = fingerprint(query_smiles)
    if query is None:
        return 0.0
    scores = []
    for reference in references:
        reference_fp = fingerprint(reference)
        if reference_fp is not None:
            scores.append(float(DataStructs.TanimotoSimilarity(query, reference_fp)))
    return max(scores, default=0.0)
