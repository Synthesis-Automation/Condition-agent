"""
Coupling reaction confirmation by attachment site analysis.

Validates cross-coupling products by connecting reactive sites and matching products.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.util.rdkit_helpers import rdkit_available
from chemtools.util.smarts_cache import compile_smarts
from chemtools.util.rdkit_helpers import parse_smiles


# SMARTS patterns for common coupling partners
_ELECTROPHILE_SMARTS = (
    "[c:1]-[F,Cl,Br,I:2]",
    "[C;X2,X3:1]-[F,Cl,Br,I:2]",
    "[c:1]-[O:2]S(=O)(=O)[#6,F]",
    "[C;X2,X3:1]-[O:2]S(=O)(=O)[#6,F]",
)
_BORON_SMARTS = (
    "[c:1]-[B:2](O)O",
    "[C;X2,X3:1]-[B:2](O)O",
    "[c:1]-[B:2](O[#6])O[#6]",
    "[C;X2,X3:1]-[B:2](O[#6])O[#6]",
    "[c:1]-[B:2]1OC(C)(C)C(C)(C)O1",
    "[C;X2,X3:1]-[B:2]1OC(C)(C)C(C)(C)O1",
    "[c:1]-[B-:2](F)(F)F",
    "[C;X2,X3:1]-[B-:2](F)(F)F",
)
_ORGANOZINC_SMARTS = (
    "[c:1]-[Zn:2]",
    "[C;X2,X3:1]-[Zn:2]",
)
_ORGANOTIN_SMARTS = (
    "[c:1]-[Sn:2]",
    "[C;X2,X3:1]-[Sn:2]",
)
_ORGANOMAGNESIUM_SMARTS = (
    "[c:1]-[Mg:2]",
    "[C;X2,X3:1]-[Mg:2]",
)
_ORGANOSILICON_SMARTS = (
    "[c:1]-[Si:2]",
    "[C;X2,X3:1]-[Si:2]",
)
_AMINE_NH_SMARTS = (
    "[N:1]",
    "[n:1]",
)
_ALCOHOL_OH_SMARTS = ("[O:1;H1;!$(O-[C,S,P]=O)]",)
_THIOL_SH_SMARTS = ("[S:1;H1;!$(S-[C,S,P]=O)]",)
_TERMINAL_ALKYNE_SMARTS = ("[C:1]#[CH]", "[C:1]#[CH1]")


@dataclass(frozen=True)
class CouplingConfirmationSpec:
    """Specification for coupling reaction confirmation."""
    electrophile_smarts: Tuple[str, ...]
    nucleophile_smarts: Tuple[str, ...]
    require_distinct_molecules: bool = True


# Registry of coupling reaction confirmation specs
COUPLING_CONFIRMATION_SPECS: Dict[str, CouplingConfirmationSpec] = {
    "Suzuki_miyaura": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _BORON_SMARTS),
    "Negishi": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _ORGANOZINC_SMARTS),
    "Stille": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _ORGANOTIN_SMARTS),
    "Kumada": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _ORGANOMAGNESIUM_SMARTS),
    "Hiyama": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _ORGANOSILICON_SMARTS),
    "Sonogashira": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _TERMINAL_ALKYNE_SMARTS),
    "C_N_Coupling": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _AMINE_NH_SMARTS),
    "C_O_Coupling": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _ALCOHOL_OH_SMARTS),
    "C_S_Coupling": CouplingConfirmationSpec(_ELECTROPHILE_SMARTS, _THIOL_SH_SMARTS),
}


def _compile_mapped_patterns(smarts_list: Iterable[str]) -> List[Tuple[Any, int, Optional[int]]]:
    """Compile SMARTS patterns and extract atom map positions."""
    compiled: List[Tuple[Any, int, Optional[int]]] = []
    for smarts in smarts_list:
        pattern = compile_smarts(smarts, validate=False)
        if not pattern:
            continue
        map_pos: Dict[int, int] = {}
        for idx, atom in enumerate(pattern.GetAtoms()):
            map_num = atom.GetAtomMapNum()
            if map_num:
                map_pos[map_num] = idx
        if 1 not in map_pos:
            continue
        compiled.append((pattern, map_pos[1], map_pos.get(2)))
    return compiled


def _find_attachment_sites(
    mol: Any,
    patterns: Iterable[str],
) -> List[Tuple[int, Optional[int]]]:
    """Find attachment sites using SMARTS patterns with atom mapping."""
    sites: Set[Tuple[int, Optional[int]]] = set()
    for pattern, map1, map2 in _compile_mapped_patterns(patterns):
        try:
            matches = mol.GetSubstructMatches(pattern)
        except Exception:
            continue
        for match in matches:
            if map1 >= len(match):
                continue
            leaving_idx = None
            if map2 is not None and map2 < len(match):
                leaving_idx = match[map2]
            sites.add((match[map1], leaving_idx))
    return sorted(sites)


def _collect_branch_atoms(mol: Any, start_idx: int, blocked_idx: int) -> Set[int]:
    """Collect all atoms in the branch starting from start_idx, stopping at blocked_idx."""
    stack = [start_idx]
    visited: Set[int] = set()
    while stack:
        idx = stack.pop()
        if idx in visited:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx == blocked_idx:
                continue
            if n_idx not in visited:
                stack.append(n_idx)
    return visited


def _adjust_index(original: int, removed: List[int]) -> int:
    """Adjust atom index after removing atoms."""
    return original - sum(1 for idx in removed if idx < original)


def _build_coupling_candidate(
    *,
    electrophile_mol: Any,
    electrophile_aryl_idx: int,
    electrophile_lg_idx: int,
    boron_mol: Any,
    boron_aryl_idx: int,
    boron_atom_idx: Optional[int],
) -> Optional[Any]:
    """Build a candidate coupling product by connecting attachment points."""
    from rdkit import Chem

    combined = Chem.CombineMols(electrophile_mol, boron_mol)
    rw = Chem.RWMol(combined)
    offset = electrophile_mol.GetNumAtoms()
    aryl_a = electrophile_aryl_idx
    aryl_b = boron_aryl_idx + offset
    lg_idx = electrophile_lg_idx
    b_idx = boron_atom_idx + offset if boron_atom_idx is not None else None

    if rw.GetBondBetweenAtoms(aryl_a, aryl_b) is None:
        rw.AddBond(aryl_a, aryl_b, Chem.BondType.SINGLE)

    to_remove = set()
    to_remove.update(_collect_branch_atoms(rw, lg_idx, aryl_a))
    if b_idx is not None:
        to_remove.update(_collect_branch_atoms(rw, b_idx, aryl_b))
    if aryl_a in to_remove or aryl_b in to_remove:
        return None
    removed = sorted(to_remove, reverse=True)
    for idx in removed:
        rw.RemoveAtom(idx)

    aryl_a = _adjust_index(aryl_a, removed)
    aryl_b = _adjust_index(aryl_b, removed)

    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None

    try:
        frags = Chem.GetMolFrags(mol, asMols=False)
    except Exception:
        return None
    for frag in frags:
        if aryl_a in frag and aryl_b in frag:
            return Chem.PathToSubmol(mol, frag)
    return None


def confirm_coupling_product_by_attachment(
    reactant_smiles: Iterable[str],
    product_smiles: Iterable[str],
    reaction_type: str,
) -> Tuple[bool, Optional[str]]:
    """Confirm a coupling product by connecting attachment-point atoms and matching products."""
    if not rdkit_available():
        return False, "rdkit_unavailable"

    spec = COUPLING_CONFIRMATION_SPECS.get(reaction_type)
    if spec is None:
        return False, "unsupported_reaction"

    reactant_mols: List[Any] = []
    for smi in reactant_smiles:
        mol = parse_smiles(smi)
        if mol is not None:
            reactant_mols.append(mol)

    product_mols: List[Any] = []
    for smi in product_smiles:
        mol = parse_smiles(smi)
        if mol is not None:
            product_mols.append(mol)

    if not reactant_mols or not product_mols:
        return False, "no_mols"

    electrophile_sites: List[Tuple[int, int, int]] = []
    nucleophile_sites: List[Tuple[int, int, Optional[int]]] = []
    for idx, mol in enumerate(reactant_mols):
        for attach_idx, lg_idx in _find_attachment_sites(mol, spec.electrophile_smarts):
            if lg_idx is None:
                continue
            electrophile_sites.append((idx, attach_idx, lg_idx))
        for attach_idx, lg_idx in _find_attachment_sites(mol, spec.nucleophile_smarts):
            nucleophile_sites.append((idx, attach_idx, lg_idx))

    if len(electrophile_sites) != 1 or len(nucleophile_sites) != 1:
        return False, "ambiguous_sites"

    e_idx, e_attach, e_lg = electrophile_sites[0]
    n_idx, n_attach, n_lg = nucleophile_sites[0]
    if spec.require_distinct_molecules and e_idx == n_idx:
        return False, "same_molecule"

    candidate = _build_coupling_candidate(
        electrophile_mol=reactant_mols[e_idx],
        electrophile_aryl_idx=e_attach,
        electrophile_lg_idx=e_lg,
        boron_mol=reactant_mols[n_idx],
        boron_aryl_idx=n_attach,
        boron_atom_idx=n_lg,
    )
    if candidate is None:
        return False, "candidate_build_failed"

    for product in product_mols:
        try:
            if product.HasSubstructMatch(candidate):
                return True, "substructure_match"
        except Exception:
            continue
    return False, "no_substructure_match"


def confirm_suzuki_product_by_attachment(
    reactant_smiles: Iterable[str],
    product_smiles: Iterable[str],
) -> Tuple[bool, Optional[str]]:
    """Backward-compatible wrapper for Suzuki confirmation."""
    return confirm_coupling_product_by_attachment(
        reactant_smiles,
        product_smiles,
        "Suzuki_miyaura",
    )
