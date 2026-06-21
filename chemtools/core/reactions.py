from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional


@dataclass(frozen=True)
class ReactionComponent:
    """Normalized reaction-side component."""

    input: str = ""
    fragments: List[str] = field(default_factory=list)
    largest_smiles: str = ""
    smiles_norm: Optional[str] = None
    error: Optional[str] = None

    @classmethod
    def from_payload(cls, payload: Any) -> "ReactionComponent":
        if not isinstance(payload, dict):
            return cls(input=str(payload or ""))
        fragments_raw = payload.get("fragments") or []
        fragments = [str(item) for item in fragments_raw if str(item)]
        return cls(
            input=str(payload.get("input") or ""),
            fragments=fragments,
            largest_smiles=str(payload.get("largest_smiles") or ""),
            smiles_norm=(
                str(payload.get("smiles_norm"))
                if payload.get("smiles_norm") is not None
                else None
            ),
            error=str(payload.get("error")) if payload.get("error") is not None else None,
        )

    def to_payload(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "input": self.input,
            "fragments": list(self.fragments),
            "largest_smiles": self.largest_smiles,
            "smiles_norm": self.smiles_norm,
        }
        if self.error is not None:
            payload["error"] = self.error
        return payload

    @property
    def preferred_smiles(self) -> str:
        return self.smiles_norm or self.largest_smiles or self.input or ""


@dataclass(frozen=True)
class ReactionRecord:
    """Structured representation of a normalized reaction."""

    input: str = ""
    reactants: List[ReactionComponent] = field(default_factory=list)
    agents: List[ReactionComponent] = field(default_factory=list)
    products: List[ReactionComponent] = field(default_factory=list)

    @staticmethod
    def _parse_side(side_payload: Iterable[Any]) -> List[ReactionComponent]:
        return [ReactionComponent.from_payload(item) for item in side_payload or []]

    @classmethod
    def from_payload(cls, payload: Dict[str, Any]) -> "ReactionRecord":
        payload = payload or {}
        return cls(
            input=str(payload.get("input") or ""),
            reactants=cls._parse_side(payload.get("reactants") or []),
            agents=cls._parse_side(payload.get("agents") or []),
            products=cls._parse_side(payload.get("products") or []),
        )

    @classmethod
    def from_component_payloads(
        cls,
        *,
        input_smiles: str,
        reactants: List[Dict[str, Any]],
        agents: List[Dict[str, Any]],
        products: List[Dict[str, Any]],
    ) -> "ReactionRecord":
        return cls(
            input=str(input_smiles or ""),
            reactants=cls._parse_side(reactants),
            agents=cls._parse_side(agents),
            products=cls._parse_side(products),
        )

    @staticmethod
    def _join_side(side: Iterable[ReactionComponent]) -> str:
        values = [component.preferred_smiles for component in side if component.preferred_smiles]
        return ".".join(values)

    @property
    def normalized(self) -> str:
        return ">".join(
            [
                self._join_side(self.reactants),
                self._join_side(self.agents),
                self._join_side(self.products),
            ]
        )

    @property
    def errors(self) -> List[str]:
        items = self.reactants + self.agents + self.products
        return [component.input for component in items if component.error]

    @property
    def reactant_smiles(self) -> List[str]:
        return [component.preferred_smiles for component in self.reactants if component.preferred_smiles]

    @property
    def product_smiles(self) -> List[str]:
        return [component.preferred_smiles for component in self.products if component.preferred_smiles]

    @property
    def agent_payloads(self) -> List[Dict[str, Any]]:
        return [component.to_payload() for component in self.agents]

    def to_payload(self) -> Dict[str, Any]:
        return {
            "input": self.input,
            "reactants": [component.to_payload() for component in self.reactants],
            "agents": [component.to_payload() for component in self.agents],
            "products": [component.to_payload() for component in self.products],
            "normalized": self.normalized,
            "errors": self.errors,
        }

