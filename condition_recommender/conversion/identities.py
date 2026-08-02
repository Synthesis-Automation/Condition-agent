"""Canonical identities for reactions, recipes, and source observations."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Iterable, Optional, Tuple

from reactive_taxonomy.chemistry.rdkit_utils import parse_smiles

from .input_schema import RawReactionRecord


def _digest(prefix: str, payload: Any) -> str:
    encoded = json.dumps(
        payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return f"{prefix}:" + hashlib.sha256(encoded).hexdigest()


@lru_cache(maxsize=65536)
def _canonical_component(smiles: str) -> Optional[str]:
    from rdkit import Chem

    mol = parse_smiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    try:
        return str(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True))
    except Exception:
        return None


def _canonical_side(side: str) -> Optional[Tuple[str, ...]]:
    values = []
    for token in (item.strip() for item in side.split(".")):
        if not token:
            continue
        canonical = _canonical_component(token)
        if canonical is None:
            return None
        values.append(canonical)
    return tuple(sorted(values))


@dataclass(frozen=True)
class CanonicalReactionIdentity:
    reaction_id: str
    canonical_reaction_smiles: str
    reactants: Tuple[str, ...]
    agents: Tuple[str, ...]
    products: Tuple[str, ...]
    schema_version: str = "1.0"


def canonical_reaction_identity(
    reaction_smiles: str,
) -> Optional[CanonicalReactionIdentity]:
    """Canonicalize unmapped components with order-invariant reaction sides."""
    text = str(reaction_smiles or "").strip()
    if ">>" in text:
        left, right = text.split(">>", 1)
        middle = ""
    else:
        parts = text.split(">")
        if len(parts) != 3:
            return None
        left, middle, right = parts
    reactants = _canonical_side(left)
    agents = _canonical_side(middle)
    products = _canonical_side(right)
    if not reactants or not products or agents is None:
        return None
    canonical = f"{'.'.join(reactants)}>{'.'.join(agents)}>{'.'.join(products)}"
    return CanonicalReactionIdentity(
        reaction_id=_digest("CRX1", {"reactants": reactants, "products": products}),
        canonical_reaction_smiles=canonical,
        reactants=reactants,
        agents=agents,
        products=products,
    )


def raw_recipe_id(record: RawReactionRecord) -> str:
    """Identify source condition groups without discarding unresolved values."""
    payload = {
        "catalyst": tuple(sorted(set(record.catalyst_cas))),
        "reagent": tuple(sorted(set(record.reagent_cas))),
        "solvent": tuple(sorted(set(record.solvent_cas))),
    }
    if (
        record.condition_component_inputs
        or record.condition_process_stages
        or record.condition_declared_absences
    ):
        payload["components"] = tuple(
            sorted(
                (
                    item.raw_identifier,
                    item.source_field,
                    item.identifier_type,
                    item.source_role_hint,
                    item.amount,
                    item.amount_unit,
                )
                for item in record.condition_component_inputs
            )
        )
        payload["stages"] = tuple(
            (
                item.stage_index,
                item.temperature_c,
                item.time_h,
                item.atmosphere,
            )
            for item in record.condition_process_stages
        )
        payload["declared_absences"] = record.condition_declared_absences
    return _digest("RAWCOND1", payload)


def observation_id(record: RawReactionRecord) -> str:
    """Identify one source observation independently of source reaction IDs."""
    return _digest(
        "OBS1",
        {
            "dataset": record.source_dataset,
            "path": record.source_path,
            "row": record.source_row_number,
        },
    )


def normalized_identifier_tokens(values: Iterable[str]) -> Tuple[str, ...]:
    return tuple(sorted({str(value).strip() for value in values if str(value).strip()}))


__all__ = [
    "CanonicalReactionIdentity",
    "canonical_reaction_identity",
    "normalized_identifier_tokens",
    "observation_id",
    "raw_recipe_id",
]
