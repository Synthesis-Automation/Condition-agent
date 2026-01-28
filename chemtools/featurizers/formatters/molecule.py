"""
Molecule-level feature formatting and bundling.

Handles molecule featurization, RDKit property calculation, and SNAr feasibility.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from chemtools.util import rdkit_helpers
from ..molecule import featurize_molecule as _featurize_molecule
from ..analysis.feasibility import analyze_molecule_snar_feasibility
from .simplified import build_core_molecule, build_extended_molecule


def to_bool(value: Any, *, default: bool) -> bool:
    """Convert any value to boolean with fallback."""
    if value is None:
        return default
    return bool(value)


def calculate_rdkit_properties(smiles: str) -> Dict[str, Any]:
    """Calculate standard RDKit molecular descriptors."""
    mol = rdkit_helpers.parse_smiles(smiles)
    if mol is None:
        return {}
    try:
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    except Exception:
        return {}
    return {
        "molecular_weight": float(Descriptors.MolWt(mol)),
        "logP": float(Descriptors.MolLogP(mol)),
        "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
        "HBA": float(Lipinski.NumHAcceptors(mol)),
        "HBD": float(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": float(Lipinski.NumRotatableBonds(mol)),
        "fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_ring_count": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "heavy_atom_count": int(Descriptors.HeavyAtomCount(mol)),
    }


def build_molecule_bundle(
    smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Build a complete molecule feature bundle.
    
    Args:
        smiles: Molecule SMILES string
        registry_paths: Custom taxonomy paths
        options: Featurization options
        
    Returns:
        Bundle with motifs, properties, steric/electronic features
    """
    options = options or {}
    include_rdkit = to_bool(options.get("include_rdkit"), default=True)
    
    result = _featurize_molecule(smiles, registry_paths=registry_paths, options=options)
    
    if include_rdkit:
        rdkit_props = calculate_rdkit_properties(smiles)
        if rdkit_props:
            result["rdkit_properties"] = rdkit_props
    
    return result


def featurize_molecule(
    smiles: str,
    *,
    registry_paths: Optional[Dict[str, str | Path]] = None,
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Return a molecule feature bundle (core or extended format).
    
    Args:
        smiles: Molecule SMILES string
        registry_paths: Custom taxonomy paths
        options: Featurization options
            - detailed (bool): If True, return extended format with all analysis.
                             If False (default), return simplified core format.
        
    Returns:
        Core bundle (6 fields) or extended bundle (7 fields with extended section)
    """
    options = options or {}
    use_detailed = to_bool(options.get("detailed"), default=False)
    
    # Build full molecule bundle first
    molecule = build_molecule_bundle(
        smiles,
        registry_paths=registry_paths,
        options=options,
    )
    
    # Add SNAr feasibility check
    snar_feasibility = analyze_molecule_snar_feasibility(molecule)
    if snar_feasibility:
        molecule["snar_feasibility"] = snar_feasibility
    
    # Return simplified format (core or extended)
    if use_detailed:
        result = build_extended_molecule(molecule)
    else:
        result = build_core_molecule(molecule)
    
    # Add metadata
    result["kind"] = "molecule"
    result["schema_version"] = "v2"
    
    meta = {
        "rdkit_available": rdkit_helpers.rdkit_available(),
    }
    if molecule.get("meta", {}).get("error"):
        meta["errors"] = [molecule["meta"]["error"]]
    if meta.get("errors") or not meta.get("rdkit_available", True):
        result["meta"] = meta
    
    return result
