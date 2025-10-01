"""Reaction normalisation wrapper for Condition MCP."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field

from chemtools.smiles import normalize_reaction as _normalize_reaction

from .base import SchemaStamped, flatten_smiles_block, validate_payload


class NormalizeReactionInput(BaseModel):
    """Input payload for the ``normalize_reaction`` tool."""

    smiles_rxn: str = Field(..., description="Reaction SMILES string (reactants>agents>products)")
    include_agents: bool = Field(
        default=True,
        description="When true, include the agent list in the response payload.",
    )


class NormalizeReactionOutput(SchemaStamped):
    """Canonicalised view of a reaction SMILES string."""

    input: str
    normalized: str
    reactants: List[str]
    products: List[str]
    agents: Optional[List[str]] = Field(default=None)
    errors: List[str] = Field(default_factory=list)


def normalize_reaction(data: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize a reaction SMILES string and return a serialisable payload."""

    payload = validate_payload(NormalizeReactionInput, data)
    raw = _normalize_reaction(payload.smiles_rxn)

    reactants = flatten_smiles_block(raw.get("reactants", []))
    products = flatten_smiles_block(raw.get("products", []))
    agents_block = flatten_smiles_block(raw.get("agents", [])) if payload.include_agents else None

    output = NormalizeReactionOutput(
        input=raw.get("input", payload.smiles_rxn),
        normalized=raw.get("normalized", payload.smiles_rxn),
        reactants=reactants,
        products=products,
        agents=agents_block,
        errors=list(raw.get("errors") or []),
    )
    return output.model_dump()
