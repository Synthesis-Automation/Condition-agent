"""Scientific workspaces for external agents, without an internal LLM controller."""

from .core.operation_contracts import OperationDefinition, OperationProvider
from .core.store import InvestigationEvent, InvestigationStore
from .workspace import ScientificWorkspace

__all__ = [
    "InvestigationEvent", "InvestigationStore", "OperationDefinition", "OperationProvider",
    "ScientificWorkspace",
]
