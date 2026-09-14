"""Graph-only identities for reaction sites and retained precursor fragments.

These descriptors never create precursor molecules or infer correspondence.
In particular, a retained aromatic fragment need not be a sanitizable molecule.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from functools import lru_cache
from typing import Any, Mapping

from rdkit import Chem

from .reaction_edits import normalize_mapped_edits
from .reaction_parser import parse_reaction_smiles
from .shared_core_graph import atom_label, canonical_graph


IDENTITY_GRAPH_VERSION = "reaction_identity_graphs.v1"
_EDIT_TYPES = frozenset({"formed", "broken", "order_changed", "hydrogen_change"})


def _digest(namespace: str, graph: str) -> str:
    payload = json.dumps([IDENTITY_GRAPH_VERSION, graph], separators=(",", ":"))
    return namespace + ":" + hashlib.sha256(payload.encode()).hexdigest()


def product_edit_site_graph(
    product: Chem.Mol,
    edits: tuple[Mapping[str, Any], ...],
) -> str:
    """Identify an anchored product site, including cleavage and H changes.

    The complete product graph preserves joint symmetry. Precursor-only
    endpoints are typed as external ports, without naming their leaving group.
    This is a site descriptor; transformation identity remains a separate key.
    """

    map_to_index = {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in product.GetAtoms()
        if atom.GetAtomMapNum() > 0
    }
    if len(map_to_index) != product.GetNumAtoms():
        return ""
    labels = [["product", atom_label(atom)] for atom in product.GetAtoms()]
    edges = [
        (
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            ["product_bond", str(bond.GetBondType()), str(bond.GetStereo())],
        )
        for bond in product.GetBonds()
    ]
    count = 0
    for edit in edits:
        kind = edit.get("edit_type")
        if kind not in _EDIT_TYPES:
            continue
        endpoints = []
        external_count = 0
        for field in ("atom_1", "atom_2"):
            atom = edit.get(field)
            if atom is None:
                continue
            index = map_to_index.get(int(atom.get("atom_map_number") or 0))
            if index is None:
                external_count += 1
            else:
                endpoints.append(index)
        if not endpoints:
            continue
        event = len(labels)
        labels.append(
            ["edit", kind, edit.get("old_order"), edit.get("new_order"), external_count]
        )
        edges.extend((index, event, "anchor") for index in endpoints)
        count += 1
    return _digest("SITE2", canonical_graph(labels, edges)) if count else ""


def retained_precursor_graph(
    reactants: Chem.Mol,
    product: Chem.Mol,
    operator_signature: str,
) -> str:
    """Describe retained fragments without sanitizing or hydrogen-capping them.

    Precursor atom states and bonds are observations. Component membership is
    retained explicitly even when deleting departure atoms disconnects a graph.
    The returned SYN2 is an identity, never an executable precursor SMILES.
    """

    product_maps = {atom.GetAtomMapNum() for atom in product.GetAtoms()}
    if 0 in product_maps or len(product_maps) != product.GetNumAtoms():
        return ""
    values = []
    for component in Chem.GetMolFrags(reactants):
        retained = [
            index
            for index in component
            if reactants.GetAtomWithIdx(index).GetAtomMapNum() in product_maps
        ]
        if not retained:
            continue
        lookup = {index: ordinal for ordinal, index in enumerate(retained)}
        labels = [atom_label(reactants.GetAtomWithIdx(index)) for index in retained]
        edges = [
            (
                lookup[bond.GetBeginAtomIdx()],
                lookup[bond.GetEndAtomIdx()],
                [str(bond.GetBondType()), str(bond.GetStereo())],
            )
            for bond in reactants.GetBonds()
            if bond.GetBeginAtomIdx() in lookup and bond.GetEndAtomIdx() in lookup
        ]
        values.append(canonical_graph(labels, edges))
    if not values:
        return ""
    return _digest(
        "SYN2", json.dumps([operator_signature, sorted(values)], separators=(",", ":"))
    )


@lru_cache(maxsize=20_000)
def provisional_product_site(reaction_smiles: str) -> str:
    """Return a supplied-map site for scheduling, without reaction admission."""

    parsed = parse_reaction_smiles(
        reaction_smiles, include_molecular_interpretation=False
    )
    if not parsed.valid or len(parsed.products) != 1:
        return ""
    normalized = normalize_mapped_edits(parsed.reactants, parsed.products)
    if not normalized.valid:
        return ""
    product = Chem.MolFromSmiles(parsed.products[0].input_smiles)
    if product is None:
        return ""
    return product_edit_site_graph(
        product, tuple(asdict(edit) for edit in normalized.edits)
    )
