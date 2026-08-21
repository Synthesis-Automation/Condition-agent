"""Stateless application facade for CLI and API assistance sessions."""

from __future__ import annotations

from typing import Any, Dict, Mapping

from .contracts import AssistanceRequest, session_from_dict
from .controller import AssistanceController
from .rendering import render_assistance


class AssistanceApplicationService:
    """Run or resume a session without persisting reaction history by default."""

    def __init__(self, controller: AssistanceController) -> None:
        self._controller = controller

    def start(
        self,
        *,
        objective: str,
        mode: str,
        structure_input: str,
        provider_settings: Mapping[str, Any] | None = None,
    ) -> Dict[str, Any]:
        """Start one explicitly advisory session."""

        request = AssistanceRequest(
            objective=objective,
            mode=mode,  # type: ignore[arg-type]
            structure_input=structure_input,
            provider_settings=dict(provider_settings or {}),
        )
        run = self._controller.run(request)
        payload = run.to_dict()
        payload["answer"] = render_assistance(run.state)
        payload["rollout_state"] = "advisory"
        return payload

    def confirm_condition_constraint(
        self,
        *,
        state: Mapping[str, Any],
        raw_value: Any,
    ) -> Dict[str, Any]:
        """Resume exactly one pending condition clarification."""

        restored = session_from_dict(state)
        run = self._controller.resume_with_condition_constraint(restored, raw_value)
        payload = run.to_dict()
        payload["answer"] = render_assistance(run.state)
        payload["rollout_state"] = "advisory"
        return payload
