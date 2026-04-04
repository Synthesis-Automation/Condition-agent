"""
Shared helper utilities for all ChemCoworker tool wrappers.

Follows the same pattern as reaction_agent/reasoning_tools.py and the
repository's LLM-facing tool wrappers for consistency.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

# Matches whitespace around the >> separator in reaction SMILES
_RXN_SEP_RE = re.compile(r'\s*>>\s*')


def _clean_rxn_smiles(smiles: str) -> str:
    """Remove whitespace around >> in a reaction SMILES string.

    Valid SMILES strings never contain spaces; this strips any stray
    whitespace that may have been introduced by LLM output or copy-paste
    (e.g. "Brc1ccccc1 >> OB(O)c1ccccc1" → "Brc1ccccc1>>OB(O)c1ccccc1").
    """
    if not smiles:
        return smiles
    return _RXN_SEP_RE.sub('>>', smiles.strip())


def _validate_reaction_smiles(
    reaction_smiles: str,
    *,
    require_product: bool = True,
) -> Tuple[str, Optional[str]]:
    """Validate a reaction SMILES string at a basic chemistry-sanity level.

    Validation rules:
    - Must contain a reaction arrow `>>`
    - Reactant side must be non-empty and parseable
    - Product side must be non-empty and parseable when `require_product=True`
    - Dot-separated components on each side must each parse as valid SMILES

    Returns:
        (cleaned_reaction_smiles, error_message_or_None)
    """
    rxn = _clean_rxn_smiles(str(reaction_smiles or ""))
    if not rxn:
        return rxn, "reaction_smiles cannot be empty"
    if ">>" not in rxn:
        return rxn, "reaction_smiles must be in 'reactants>>products' format"

    left, right = rxn.split(">>", 1)
    left = left.strip()
    right = right.strip()
    if not left:
        return rxn, "reaction_smiles is missing reactants before '>>'"
    if require_product and not right:
        return rxn, "reaction_smiles is missing product after '>>'"

    try:
        from rdkit import Chem, rdBase
    except Exception:
        # RDKit is required in this project; if import fails, keep structural checks only.
        return rxn, None

    def _side_ok(side: str) -> bool:
        parts = [p.strip() for p in side.split(".") if p.strip()]
        if not parts:
            return False
        for part in parts:
            with rdBase.BlockLogs():
                mol = Chem.MolFromSmiles(part)
            if mol is None:
                return False
        return True

    if not _side_ok(left):
        return rxn, "reaction_smiles has invalid reactant component(s)"
    if right and not _side_ok(right):
        return rxn, "reaction_smiles has invalid product component(s)"

    return rxn, None


def _to_jsonable(value: Any) -> Any:
    """Recursively convert to JSON-serializable types."""
    from dataclasses import asdict, is_dataclass
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return {k: _to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        return _to_jsonable(to_dict())
    to_list = getattr(value, "tolist", None)
    if callable(to_list):
        return _to_jsonable(to_list())
    return value


def _success(data: Any) -> Dict[str, Any]:
    """Wrap a successful result."""
    payload: Dict[str, Any] = {"success": True}
    if isinstance(data, dict):
        payload.update(_to_jsonable(data))
    else:
        payload["data"] = _to_jsonable(data)
    return payload


def _error(message: str, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Wrap an error result."""
    payload: Dict[str, Any] = {"success": False, "error": message}
    if extra:
        payload.update(_to_jsonable(extra))
    return payload
