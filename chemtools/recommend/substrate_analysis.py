"""
Generic substrate analysis utilities.

Functions for analyzing molecular properties and functional groups.
"""

from __future__ import annotations

from typing import Optional

from ..util.rdkit_helpers import rdkit_available, parse_smiles


def has_free_alcohol(smiles: Optional[str]) -> bool:
    """
    Check if molecule has free alcohol group (not part of carboxylic acid).
    
    Args:
        smiles: SMILES string
        
    Returns:
        True if free alcohol found
    """
    if not smiles:
        return False
    
    if rdkit_available():
        try:
            from rdkit import Chem  # type: ignore
        except Exception:
            Chem = None  # type: ignore
        
        if Chem is not None:
            mol = parse_smiles(smiles)
            if mol is not None:
                try:
                    alcohol = Chem.MolFromSmarts("[OX2H]")
                    carboxy_acid = Chem.MolFromSmarts("C(=O)[OX2H]")
                    acid_matches = {match[2] for match in mol.GetSubstructMatches(carboxy_acid)} if carboxy_acid else set()
                    
                    for match in mol.GetSubstructMatches(alcohol):
                        oxygen_idx = match[0]
                        if oxygen_idx in acid_matches:
                            continue
                        
                        atom = mol.GetAtomWithIdx(oxygen_idx)
                        is_carboxylic = False
                        for nbr in atom.GetNeighbors():
                            bond = mol.GetBondBetweenAtoms(oxygen_idx, nbr.GetIdx())
                            if nbr.GetAtomicNum() == 6 and bond is not None:
                                for b in nbr.GetBonds():
                                    if b.GetOtherAtom(nbr).GetIdx() == oxygen_idx:
                                        continue
                                    if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(nbr).GetAtomicNum() == 8:
                                        is_carboxylic = True
                                        break
                                if is_carboxylic:
                                    break
                        if is_carboxylic:
                            continue
                        return True
                except Exception:
                    pass
    
    # Fallback: simple string matching
    s_lower = smiles.lower()
    if "oh" in s_lower or "[oh]" in s_lower:
        if "c(=o)o" in s_lower:
            return "oo" in s_lower or s_lower.count("oh") > 1
        return True
    return False


def has_phenol(smiles: Optional[str]) -> bool:
    """
    Check if molecule has phenol group.
    
    Args:
        smiles: SMILES string
        
    Returns:
        True if phenol found
    """
    if not smiles:
        return False
    
    if rdkit_available():
        try:
            from rdkit import Chem  # type: ignore
            mol = parse_smiles(smiles)
            if mol is not None:
                phenol_pattern = Chem.MolFromSmarts("c[OX2H]")
                if phenol_pattern and mol.HasSubstructMatch(phenol_pattern):
                    return True
        except Exception:
            pass
    
    # Fallback: simple pattern
    s_lower = smiles.lower()
    return "coh" in s_lower or "c[oh]" in s_lower


def has_sulfonamide(smiles: Optional[str]) -> bool:
    """
    Check if molecule has sulfonamide group.
    
    Args:
        smiles: SMILES string
        
    Returns:
        True if sulfonamide found
    """
    if not smiles:
        return False
    
    s_lower = smiles.lower()
    return "s(=o)(=o)n" in s_lower or "so2n" in s_lower


def has_hydroxylamine(smiles: Optional[str]) -> bool:
    """
    Check if molecule has hydroxylamine group.
    
    Args:
        smiles: SMILES string
        
    Returns:
        True if hydroxylamine found
    """
    if not smiles:
        return False
    
    if rdkit_available():
        try:
            from rdkit import Chem  # type: ignore
            mol = parse_smiles(smiles)
            if mol is not None:
                pattern = Chem.MolFromSmarts("[NX3][OX2H]")
                if pattern and mol.HasSubstructMatch(pattern):
                    return True
        except Exception:
            pass
    
    s_lower = smiles.lower()
    return "no" in s_lower or "n-o" in s_lower


def detect_functional_groups(smiles: Optional[str]) -> dict[str, bool]:
    """
    Detect common functional groups in molecule.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dict of functional group flags
    """
    return {
        "free_alcohol": has_free_alcohol(smiles),
        "phenol": has_phenol(smiles),
        "sulfonamide": has_sulfonamide(smiles),
        "hydroxylamine": has_hydroxylamine(smiles),
    }
