from __future__ import annotations

"""Essential-core normalization for SchemeConditionDB matching."""

from typing import List, Sequence, Set, Tuple

from rdkit import Chem

try:  # rdMolStandardize is optional in some RDKit builds
    from rdkit.Chem import rdMolStandardize
except ImportError:  # pragma: no cover - fallback path when RDKit lacks standardize module
    rdMolStandardize = None  # type: ignore[assignment]

from .types import AnchorHit, EssentialCoreNormalizationResult

# SMARTS anchors that define the core functional handles we care about.
ANCHOR_SMARTS: Tuple[str, ...] = (
    "[c]-[Cl,Br,I]",
    "[c]-OS(=O)(=O)[C,F]",
    "[NH2]-[c]",
    "[NH]([CX4])[CX4]",
    "[N;H1]C(=O)O",
)
ANCHOR_PATTERNS: Tuple[Tuple[str, Chem.Mol], ...] = tuple(
    (smarts, Chem.MolFromSmarts(smarts)) for smarts in ANCHOR_SMARTS
)

# Atomic numbers regarded as inorganic counterions for quick stripping.
_COUNTERION_ATOMIC_NUMS: Set[int] = {
    1, 3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 32, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
}


def _is_counterion_only(mol: Chem.Mol) -> bool:
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not heavy_atoms:
        return True
    return all(atom.GetAtomicNum() in _COUNTERION_ATOMIC_NUMS for atom in heavy_atoms)


def _sanitize_mol(mol: Chem.Mol) -> Chem.Mol:
    Chem.SanitizeMol(mol)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
        Chem.SetAromaticity(mol)
    except Exception:
        # Fall back to sanitized structure if kekulization fails
        pass
    return mol


def essential_core_normalize(reactant_smiles_list: Sequence[str]) -> EssentialCoreNormalizationResult:
    """Normalize reactants ahead of SMARTS-based scheme matching."""

    if not reactant_smiles_list:
        raise ValueError("No reactant SMILES provided for normalization")

    uncharger = rdMolStandardize.Uncharger() if rdMolStandardize is not None else None
    sanitized_mols: List[Chem.Mol] = []
    masked_mols: List[Chem.Mol] = []
    masked_sources: List[int] = []
    anchor_hits: List[AnchorHit] = []
    smarts_bag: Set[str] = set()
    kept_smiles: List[str] = []
    dropped_smiles: List[str] = []

    for smiles in reactant_smiles_list:
        if not smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid reactant SMILES: {smiles!r}")
        if _is_counterion_only(mol):
            dropped_smiles.append(smiles)
            continue

        if uncharger is not None:
            mol = uncharger.uncharge(mol)
        mol = _sanitize_mol(mol)
        kept_index = len(sanitized_mols)
        sanitized_mols.append(mol)
        masked_mols.append(mol)
        masked_sources.append(kept_index)
        kept_smiles.append(smiles)

        for anchor_smarts, anchor_mol in ANCHOR_PATTERNS:
            if anchor_mol is None:
                continue
            matches = mol.GetSubstructMatches(anchor_mol)
            if not matches:
                continue
            smarts_bag.add(anchor_smarts)
            for match in matches:
                anchor_hits.append(
                    AnchorHit(
                        anchor_smarts=anchor_smarts,
                        reactant_index=kept_index,
                        atom_indices=tuple(match),
                    )
                )

    if not sanitized_mols:
        raise ValueError("All reactants were dropped during normalization")

    return EssentialCoreNormalizationResult(
        reactant_smiles=tuple(reactant_smiles_list),
        kept_reactants=tuple(kept_smiles),
        dropped_reactants=tuple(dropped_smiles),
        sanitized_mols=tuple(sanitized_mols),
        masked_mols=tuple(masked_mols),
        masked_source_indices=tuple(masked_sources),
        anchor_hits=tuple(anchor_hits),
        smarts_bag=tuple(sorted(smarts_bag)),
    )
