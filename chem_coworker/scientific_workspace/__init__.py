"""Scientific workspaces for external agents, without an internal LLM controller."""

from .store import InvestigationEvent, InvestigationStore
from .operation_contracts import OperationDefinition, OperationProvider
from .workspace import ScientificWorkspace

__all__ = [
    "InvestigationEvent", "InvestigationStore", "OperationDefinition", "OperationProvider",
    "ScientificWorkspace",
]
