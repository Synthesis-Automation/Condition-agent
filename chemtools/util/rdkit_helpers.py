from typing import Optional
def canonical_smiles(smiles: str) -> Optional[str]:
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return None
    m = Chem.MolFromSmiles(smiles)
    if not m: return None
    return Chem.MolToSmiles(m, canonical=True)
