"""Stable execution contracts for independently registered scientific capabilities."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Protocol


def _identity(result: Any) -> Any:
    return result


@dataclass(frozen=True)
class OperationDefinition:
    """Declare recording policy without putting tool-specific rules in the core.

    Implementations and replay projections are registered Python code. Dataset
    requirements describe discovery prerequisites; tools still validate optional
    dependencies and scientific input themselves when invoked. A declared execution
    status field must contain completed, error, timed_out or cancelled; scientific
    uncertainty belongs in separate domain fields.
    """

    name: str
    contract_version: str = "1"
    required_artifacts: tuple[str, ...] = ()
    evidence_arguments: tuple[str, ...] = ()
    usage_policy: str = ""
    execution_status_field: str | None = None
    replay_comparison: str = "full_result"
    replay_projection: Callable[[Any], Any] = field(default=_identity, repr=False)

    def __post_init__(self) -> None:
        """Reject invalid declarations before a capability can be used."""
        if not self.name.isidentifier() or self.name.startswith("_"):
            raise ValueError("Operation name must be a public Python identifier")
        if not self.contract_version.strip():
            raise ValueError("Operation contract_version must be nonempty")
        for name in ("required_artifacts", "evidence_arguments"):
            values = getattr(self, name)
            if not isinstance(values, tuple) or any(not isinstance(value, str) or not value for value in values):
                raise ValueError(f"Operation {name} must be an immutable tuple of nonempty names")
        if not callable(self.replay_projection):
            raise ValueError("Replay projection must be registered Python code")

    def evidence_references(self, arguments: Mapping[str, Any]) -> tuple[str, ...]:
        """Find declared reference inputs; the store verifies their actual contents."""
        references: list[str] = []
        for name in self.evidence_arguments:
            value = arguments.get(name)
            if isinstance(value, str):
                references.append(value)
            elif isinstance(value, (list, tuple)):
                references.extend(item for item in value if isinstance(item, str))
        return tuple(dict.fromkeys(references))


class OperationProvider(Protocol):
    """The small discovery and invocation interface consumed by the workspace."""

    def definition(self, operation: str) -> OperationDefinition:
        """Return a declared operation or raise ValueError for an unknown name."""
        ...

    def catalog(self) -> list[dict[str, Any]]:
        """Describe available contracts, arguments, and capability prerequisites."""
        ...

    def invoke(self, operation: str, arguments: Mapping[str, Any]) -> Any:
        """Execute an explicitly registered capability."""
        ...
