from typing import Dict, Any, List
from .util.rdkit_helpers import (
    canonical_smiles,
    parse_smiles,
    neutralize_and_standardize,
    mol_to_canonical_smiles,
    choose_largest_organic_fragment,
    rdkit_available,
)

def _split_fragments_text(s: str) -> List[str]:
    if not s:
        return []
    return [p for p in s.split('.') if p]

def normalize(smi: str) -> Dict[str, Any]:
    s = (smi or '').strip()
    frags_text = _split_fragments_text(s)

    mol = parse_smiles(s)
    if mol is None:
        # If RDKit is unavailable, fallback to heuristics and do not error
        if not rdkit_available():
            largest = frags_text[0] if frags_text else ""
            # Basic textual cleanups for common salts/charges
            if largest:
                txt = largest
                for bad, rep in (
                    ('[Na+]', ''), ('[K+]', ''), ('[Li+]', ''), ('[Cl-]', 'Cl'), ('[Br-]', 'Br'), ('[I-]', 'I'),
                ):
                    txt = txt.replace(bad, rep)
                # Carboxylate to acid (common ordering variants)
                import re as _re
                txt = _re.sub(r"C\(=O\)\[O-\]", "C(=O)O", txt)
                txt = _re.sub(r"\[O-\]C\(=O\)", "OC(=O)", txt)
                txt = _re.sub(r"O=C\(\[O-\]\)", "O=C(O)", txt)
                # Simplistic ammonium to amine
                txt = txt.replace('[N+]', 'N')
                largest = txt
            return {
                "input": s,
                "fragments": frags_text,
                "largest_smiles": largest,
                "smiles_norm": largest,
            }
        # RDKit available but cannot parse -> invalid SMILES
        return {
            "input": s,
            "fragments": frags_text,
            "largest_smiles": frags_text[0] if frags_text else "",
            "smiles_norm": None,
            "error": "INVALID_SMILES",
        }
    # If RDKit present, standardize
    mol_std = neutralize_and_standardize(mol) if mol is not None else None
    if mol_std is None and mol is not None:
        mol_std = choose_largest_organic_fragment(mol)

    # Determine largest fragment string from RDKit if possible, else textual
    largest_smiles = mol_to_canonical_smiles(mol_std) if mol_std is not None else (frags_text[0] if frags_text else "")

    # Canonical normalized SMILES
    smiles_norm = mol_to_canonical_smiles(mol_std) if mol_std is not None else canonical_smiles(largest_smiles) or largest_smiles

    return {
        "input": s,
        "fragments": frags_text,
        "largest_smiles": largest_smiles,
        "smiles_norm": smiles_norm,
    }


# --- Reaction SMILES helpers ---

def _split_reaction_smiles(rsmi: str) -> List[str]:
    # Reaction SMILES has 3 fields: reactants > agents > products
    # Some strings may use '>>' with empty agents.
    parts = rsmi.split('>')
    if len(parts) == 2 and '>>' in rsmi:
        # empty agents
        return [parts[0], '', parts[1]]
    if len(parts) == 3:
        return parts
    # Fallback: treat all as reactants only
    return [rsmi, '', '']

def _normalize_list(segment: str) -> List[Dict[str, Any]]:
    if not segment:
        return []
    mols = [m for m in segment.split('.') if m]
    return [normalize(m) for m in mols]

def normalize_reaction(rsmi: str) -> Dict[str, Any]:
    r, a, p = _split_reaction_smiles((rsmi or '').strip())
    reactants = _normalize_list(r)
    agents = _normalize_list(a)
    products = _normalize_list(p)
    def join_side(items: List[Dict[str, Any]]) -> str:
        out = []
        for it in items:
            s = it.get('smiles_norm') or it.get('largest_smiles') or it.get('input') or ''
            if s:
                out.append(s)
        return '.'.join(out)
    normalized = '>'.join([join_side(reactants), join_side(agents), join_side(products)])
    errors = []
    for it in reactants + agents + products:
        if it.get('error'):
            errors.append(it['input'])
    return {
        "input": rsmi,
        "reactants": reactants,
        "agents": agents,
        "products": products,
        "normalized": normalized,
        "errors": errors,
    }
