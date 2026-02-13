"""Session memory for reaction agent runs."""

from __future__ import annotations

from dataclasses import asdict
from typing import Dict, List, Optional
from uuid import uuid4

from .contracts import SessionState, TraceEvent, utc_now_iso


class SessionMemory:
    """In-memory session store with serialized per-session state."""

    def __init__(self) -> None:
        self._sessions: Dict[str, SessionState] = {}

    def create_session(self, reaction_smiles: str, session_id: Optional[str] = None) -> SessionState:
        """Create and register a new session."""
        resolved_id = session_id or f"rxn-{uuid4().hex[:12]}"
        state = SessionState(session_id=resolved_id, reaction_smiles=reaction_smiles)
        self._sessions[resolved_id] = state
        return state

    def get_session(self, session_id: str) -> Optional[SessionState]:
        """Get an existing session by ID."""
        return self._sessions.get(session_id)

    def append_trace(
        self,
        state: SessionState,
        *,
        step_index: int,
        action: str,
        tool_name: Optional[str],
        status: str,
        message: str,
    ) -> None:
        """Append a trace event to a session."""
        state.trace.append(
            TraceEvent(
                step_index=step_index,
                timestamp=utc_now_iso(),
                action=action,
                tool_name=tool_name,
                status=status,
                message=message,
            )
        )

    def list_sessions(self) -> List[str]:
        """List session IDs."""
        return sorted(self._sessions)

    def export_session(self, session_id: str) -> Dict[str, object]:
        """Export a session as a dictionary."""
        state = self.get_session(session_id)
        if state is None:
            return {}
        return asdict(state)
