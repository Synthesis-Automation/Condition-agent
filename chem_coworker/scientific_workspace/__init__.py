"""Scientific workspaces for external agents, without an internal LLM controller."""

from .store import InvestigationEvent, InvestigationStore
from .workspace import ScientificWorkspace

__all__ = ["InvestigationEvent", "InvestigationStore", "ScientificWorkspace"]
