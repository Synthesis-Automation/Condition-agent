"""
Coupling reaction confirmation by attachment site analysis.

Validates cross-coupling products by connecting reactive sites and matching products.
All SMARTS patterns are loaded from reaction taxonomy (reaction_types.v4.0.json).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from chemtools.taxonomy import reaction_catalog
from chemtools.util.rdkit_helpers import rdkit_available, parse_smiles
from chemtools.util.smarts_cache import compile_smarts


@dataclass(frozen=True)
class CouplingConfirmationSpec:
    """Specification for coupling reaction confirmation."""
    electrophile_smarts: Tuple[str, ...]
    nucleophile_smarts: Tuple[str, ...]
    require_distinct_molecules: bool = True


@lru_cache(maxsize=1)
def _load_coupling_specs_from_taxonomy() -> Dict[str, CouplingConfirmationSpec]:
    """Load coupling confirmation specs from reaction taxonomy.
    
    Reads SMARTS patterns for electrophiles and nucleophiles from each
    coupling reaction's metadata in reaction_types.v4.0.json.
    
    Returns:
        Dictionary mapping reaction_type -> CouplingConfirmationSpec
    """
    catalog = reaction_catalog.load_reaction_catalog()
    reaction_types = catalog.get("reaction_types", [])
    
    specs: Dict[str, CouplingConfirmationSpec] = {}
    
    # Default electrophile patterns (used as fallback if not specified)
    default_electrophile_smarts = (
        "[c:1]-[F,Cl,Br,I:2]",
        "[C;X2,X3:1]-[F,Cl,Br,I:2]",
        "[c:1]-[O:2]S(=O)(=O)[#6,F]",
        "[C;X2,X3:1]-[O:2]S(=O)(=O)[#6,F]",
    )
    
    for rxn_data in reaction_types:
        if not isinstance(rxn_data, dict):
            continue
        
        rxn_id = rxn_data.get("id")
        if not rxn_id:
            continue
            
        metadata = rxn_data.get("metadata", {})
        coupling_config = metadata.get("coupling_confirmation")
        
        if not coupling_config or not isinstance(coupling_config, dict):
            continue
            
        # Extract patterns from taxonomy
        electrophile_patterns = coupling_config.get("electrophile_smarts")
        nucleophile_patterns = coupling_config.get("nucleophile_smarts")
        
        if not nucleophile_patterns:
            continue
            
        # Use default electrophiles if not specified
        if not electrophile_patterns:
            electrophile_patterns = default_electrophile_smarts
            
        # Convert to tuples
        if isinstance(electrophile_patterns, list):
            electrophile_patterns = tuple(electrophile_patterns)
        elif isinstance(electrophile_patterns, str):
            electrophile_patterns = (electrophile_patterns,)
            
        if isinstance(nucleophile_patterns, list):
            nucleophile_patterns = tuple(nucleophile_patterns)
        elif isinstance(nucleophile_patterns, str):
            nucleophile_patterns = (nucleophile_patterns,)
            
        require_distinct = coupling_config.get("require_distinct_molecules", True)
        
        specs[rxn_id] = CouplingConfirmationSpec(
            electrophile_smarts=electrophile_patterns,
            nucleophile_smarts=nucleophile_patterns,
            require_distinct_molecules=require_distinct,
        )
    
    return specs


def get_coupling_confirmation_specs() -> Dict[str, CouplingConfirmationSpec]:
    """Get all coupling confirmation specs from taxonomy (cached).
    
    Returns:
        Dictionary mapping reaction_type -> CouplingConfirmationSpec
    """
    return _load_coupling_specs_from_taxonomy()


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
    """Confirm a coupling product by connecting attachment-point atoms and matching products.
    
    Uses SMARTS patterns loaded from reaction taxonomy to identify electrophile
    and nucleophile attachment sites, then builds a candidate product and checks
    if it matches the actual product.
    
    Args:
        reactant_smiles: SMILES strings of reactants
        product_smiles: SMILES strings of products
        reaction_type: Reaction type identifier (e.g., "Suzuki_miyaura", "C_N_Coupling")
        
    Returns:
        Tuple of (success: bool, reason: str)
    """
    if not rdkit_available():
        return False, "rdkit_unavailable"

    # Load coupling specs from taxonomy
    specs = get_coupling_confirmation_specs()
    spec = specs.get(reaction_type)
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


# Legacy Suzuki-specific function removed - use confirm_coupling_product_by_attachment instead
