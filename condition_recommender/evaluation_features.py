"""Deterministic molecular features used only for leakage-safe evaluation."""

from __future__ import annotations

import hashlib
from typing import Any, Mapping, Tuple

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def reactive_component_indices(signature: Mapping[str, Any]) -> Tuple[int, ...]:
    """Return reactant component indices supported by partners or bond edits."""
    partner_indices = {
        int(partner["component_index"])
        for partner in signature.get("partners") or ()
        if isinstance(partner, Mapping)
        and partner.get("component_index") is not None
    }
    if partner_indices:
        return tuple(sorted(partner_indices))
    edit_indices = {
        int(atom["component_index"])
        for edit in signature.get("edits") or ()
        if isinstance(edit, Mapping)
        for atom in (edit.get("atom_1"), edit.get("atom_2"))
        if isinstance(atom, Mapping)
        and atom.get("side") == "reactant"
        and atom.get("component_index") is not None
    }
    return tuple(sorted(edit_indices))


def _reactant_text(reaction_smiles: str) -> str:
    if ">>" in reaction_smiles:
        return reaction_smiles.split(">>", maxsplit=1)[0]
    parts = reaction_smiles.split(">")
    return parts[0] if len(parts) == 3 else ""


def reaction_scaffold_tokens(
    reaction_smiles: str,
    signature: Mapping[str, Any],
) -> Tuple[str, ...]:
    """Return typed Bemis-Murcko tokens for structurally reactive components.

    Acyclic reactive partners use their canonical molecular graph as the strict
    leakage token because a Bemis-Murcko scaffold is empty for those molecules.
    Invalid components are omitted; conversion validity is handled upstream.
    """
    components = _reactant_text(reaction_smiles).split(".")
    indices = reactive_component_indices(signature)
    if not indices:
        indices = tuple(range(len(components)))
    tokens = []
    for index in indices:
        if index < 0 or index >= len(components):
            continue
        molecule = Chem.MolFromSmiles(components[index])
        if molecule is None:
            continue
        scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
        if scaffold.GetNumAtoms():
            smiles = Chem.MolToSmiles(
                scaffold,
                canonical=True,
                isomericSmiles=False,
            )
            tokens.append(f"BM:{smiles}")
        else:
            smiles = Chem.MolToSmiles(
                molecule,
                canonical=True,
                isomericSmiles=False,
            )
            tokens.append(f"ACYCLIC:{smiles}")
    return tuple(sorted(tokens))


def reaction_scaffold_key(
    reaction_smiles: str,
    signature: Mapping[str, Any],
) -> str:
    """Return a stable identity for the reactive scaffold combination."""
    tokens = reaction_scaffold_tokens(reaction_smiles, signature)
    if not tokens:
        return ""
    digest = hashlib.sha256("\0".join(tokens).encode("utf-8")).hexdigest()
    return f"SCF1:{digest}"


__all__ = [
    "reaction_scaffold_key",
    "reaction_scaffold_tokens",
    "reactive_component_indices",
]
