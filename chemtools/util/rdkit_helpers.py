from typing import Optional, List
import os

def _import_rdkit():
    # Allow tests or deployments to disable RDKit via env for speed/portability
    if os.environ.get("CHEMTOOLS_DISABLE_RDKIT", "").strip().lower() in {"1", "true", "yes", "on"}:
        return None, None
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return None, None
    # rdMolStandardize is optional across RDKit builds; degrade gracefully
    try:
        from rdkit.Chem import rdMolStandardize  # type: ignore
    except Exception:
        rdMolStandardize = None  # type: ignore
    return Chem, rdMolStandardize

def rdkit_available() -> bool:
    Chem, _ = _import_rdkit()
    return Chem is not None

def canonical_smiles(smiles: str) -> Optional[str]:
    Chem, _ = _import_rdkit()
    if Chem is None:
        return None
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return None
    return Chem.MolToSmiles(m, canonical=True)

def parse_smiles(smiles: str):
    Chem, _ = _import_rdkit()
    if Chem is None:
        return None
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None

def choose_largest_organic_fragment(mol):
    """Prefer the largest fragment containing carbon; fallback to largest by heavy atoms."""
    Chem, rdMolStandardize = _import_rdkit()
    if Chem is None:
        return None
    if mol is None:
        return None
    try:
        frags = Chem.GetMolFrags(mol, asMols=True)
    except Exception:
        return mol
    def score(m):
        atoms = list(m.GetAtoms())
        has_c = any(a.GetAtomicNum() == 6 for a in atoms)
        heavy = sum(1 for a in atoms if a.GetAtomicNum() > 1)
        return (1 if has_c else 0, heavy)
    best = None
    best_s = (-1, -1)
    for f in frags:
        s = score(f)
        if s > best_s:
            best = f
            best_s = s
    return best or mol

def neutralize_and_standardize(mol):
    """Apply RDKit standard cleanup, uncharge, and pick largest organic fragment."""
    Chem, rdMolStandardize = _import_rdkit()
    if Chem is None:
        return None
    if mol is None:
        return None
    try:
        # General cleanup (tautomer/charges standardization where applicable)
        mol = rdMolStandardize.Cleanup(mol)
    except Exception:
        pass
    try:
        # Remove common salts/solvents if present
        remover = rdMolStandardize.SaltRemover()
        mol = remover.StripMol(mol, dontRemoveEverything=True)
    except Exception:
        pass
    try:
        # Choose largest (organic) fragment
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)
    except Exception:
        # Fallback manual chooser
        mol = choose_largest_organic_fragment(mol)
    try:
        # Uncharge common anions/cations (carboxylate, ammonium, etc.)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
    except Exception:
        pass
    try:
        # Canonicalize tautomers to a deterministic representative
        te = rdMolStandardize.TautomerEnumerator()
        mol = te.Canonicalize(mol)
    except Exception:
        pass
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    return mol

def mol_to_canonical_smiles(mol) -> Optional[str]:
    Chem, _ = _import_rdkit()
    if Chem is None or mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None
