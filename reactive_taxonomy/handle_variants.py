"""Explicit hypothetical lookup variants at observed aromatic leaving handles."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from rdkit import Chem

from .chemistry.rdkit_utils import parse_smiles
from .reaction_parser import parse_reaction_smiles


@dataclass(frozen=True)
class RelatedHandleVariant:
    """A search hypothesis, never an observation or predicted product."""

    reaction_smiles: str
    query_element: str
    precedent_element: str
    component_index: int
    atom_index: int
    schema_version: str = "1.0"


def _key(atom: Mapping[str, Any]) -> tuple[Any, ...]:
    return atom.get("side"), atom.get("component_index"), atom.get("atom_index")


def aromatic_leaving_handle_variant(
    reaction_smiles: str,
    signature: Mapping[str, Any],
    *,
    query_element: str,
    precedent_element: str,
    formed_partner_elements: tuple[str, ...],
) -> RelatedHandleVariant | None:
    """Replace one observed departing Br/I atom for bounded analogue lookup.

    Require a single substitution event, an explicit broken aromatic C-X bond,
    and one formed bond at that same carbon. Retained halogens, multiple sites,
    order changes and contradictory atom references cannot generate variants.
    Products and agent components are retained verbatim.
    """
    if {query_element, precedent_element} != {"I", "Br"}:
        return None
    if signature.get("event_count") != 1:
        return None
    edits = tuple(signature.get("edits") or ())
    if any(edit.get("edit_type") not in {"formed", "broken", "hydrogen_change"}
           for edit in edits):
        return None
    formed = [edit for edit in edits if edit.get("edit_type") == "formed"]
    if len(formed) != 1 or formed[0].get("new_order") != "SINGLE":
        return None
    formed_atoms = [formed[0].get("atom_1") or {}, formed[0].get("atom_2") or {}]
    if any(atom.get("side") != "reactant" for atom in formed_atoms):
        return None
    candidates = []
    for edit in edits:
        if edit.get("edit_type") != "broken" or edit.get("old_order") != "SINGLE":
            continue
        endpoints = [edit.get("atom_1") or {}, edit.get("atom_2") or {}]
        for leaving, center in (endpoints, endpoints[::-1]):
            if (leaving.get("element") != query_element
                    or center.get("element") != "C" or not center.get("aromatic")
                    or leaving.get("side") != "reactant"
                    or center.get("side") != "reactant"
                    or leaving.get("component_index") != center.get("component_index")):
                continue
            partners = [atom for atom in formed_atoms if _key(atom) != _key(center)]
            if (not any(_key(atom) == _key(center) for atom in formed_atoms)
                    or len(partners) != 1
                    or partners[0].get("element") not in formed_partner_elements):
                continue
            candidates.append((leaving, center))
    if len(candidates) != 1:
        return None
    leaving, center = candidates[0]
    parsed = parse_reaction_smiles(reaction_smiles, include_molecular_interpretation=False)
    if not parsed.valid:
        return None
    components = {part.component_index: part for part in parsed.reactants}
    component = components.get(leaving.get("component_index"))
    molecule = parse_smiles(component.input_smiles) if component else None
    if molecule is None:
        return None
    root_index, center_index = leaving.get("atom_index"), center.get("atom_index")
    if (not isinstance(root_index, int) or not isinstance(center_index, int)
            or not 0 <= root_index < molecule.GetNumAtoms()
            or not 0 <= center_index < molecule.GetNumAtoms()):
        return None
    root = molecule.GetAtomWithIdx(root_index)
    anchor = molecule.GetAtomWithIdx(center_index)
    bond = molecule.GetBondBetweenAtoms(root_index, center_index)
    if (root.GetSymbol() != query_element or root.GetDegree() != 1
            or root.GetFormalCharge() != 0 or root.GetIsotope() != 0
            or root.GetNumRadicalElectrons() != 0
            or anchor.GetSymbol() != "C" or not anchor.GetIsAromatic()
            or bond is None or bond.GetBondType() != Chem.BondType.SINGLE):
        return None
    if any(
        reference.get("atom_map_number") is not None
        and reference["atom_map_number"] != atom.GetAtomMapNum()
        for reference, atom in ((leaving, root), (center, anchor))
    ):
        return None
    parts = reaction_smiles.split(">")
    if len(parts) != 3:
        return None
    if root.GetAtomMapNum():
        product = parse_smiles(parts[2])
        if product is None or any(
            atom.GetAtomMapNum() == root.GetAtomMapNum() for atom in product.GetAtoms()
        ):
            return None
    edited = Chem.RWMol(molecule)
    edited.GetAtomWithIdx(root_index).SetAtomicNum(
        Chem.GetPeriodicTable().GetAtomicNumber(precedent_element)
    )
    try:
        Chem.SanitizeMol(edited)
    except (ValueError, RuntimeError):
        return None
    reactants = [
        Chem.MolToSmiles(edited) if part.component_index == component.component_index
        else part.input_smiles for part in parsed.reactants
    ]
    return RelatedHandleVariant(
        reaction_smiles=f"{'.'.join(reactants)}>{parts[1]}>{parts[2]}",
        query_element=query_element,
        precedent_element=precedent_element,
        component_index=component.component_index,
        atom_index=root_index,
    )
