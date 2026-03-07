"""
Phase 3/6: WorkflowRegistry — replaces hardcoded if task_type == 'retrosynthesis' routing.
Phase 6 adds CriticStep: adversarial reviewer for retrosynthesis.

Usage:
    from chem_coworker.workflow import WORKFLOW_REGISTRY

    workflow = WORKFLOW_REGISTRY.get_for_task("retrosynthesis")
    print(workflow.name, workflow.max_iterations)
    # → retrosynthesis 10

    # In agent.run(), if critic step is configured:
    if workflow.critic_step and workflow.critic_step.enabled:
        findings, verdict = critic_agent.review(...)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core types
# ---------------------------------------------------------------------------

@dataclass
class CriticStep:
    """
    Configuration for the adversarial critic pass on a workflow.

    Fields
    ------
    enabled : bool
        If False, the critic is never called even if a CriticStep is present.
    max_findings : int
        Cap on findings returned by CriticAgent.review().
    min_severity : str
        Minimum severity level to surface: "info", "warning", or "critical".
        Findings below this level are silently dropped.
    revision_pass : bool
        If True, run a follow-up LLM call after the critic to produce a revised
        answer that addresses the identified findings.  Only fires when at least
        one finding meets the min_severity threshold.
    """
    enabled: bool = True
    max_findings: int = 5
    min_severity: str = "warning"   # mirrors critic.Severity without a circular import
    revision_pass: bool = True


@dataclass
class WorkflowDefinition:
    """
    All the settings that make one task type behave differently from another.

    Fields
    ------
    name : str
        Canonical slug — e.g. "retrosynthesis", "forward_chemistry".
    system_prompt : str
        Static SystemMessage content (can be cached by the API).
    classifier_predicate : Callable[[str], bool]
        Given a task_type string, returns True if this workflow handles it.
        First matching definition wins; the fallback definition should always return True.
    max_iterations : int
        Hard cap on native tool-loop iterations (default: 8).
    critic_step : optional
        CriticStep dataclass set in Phase 6; None until then.
    llm_visible_tools : optional
        Curated subset of registered tools exposed to the LLM for this workflow.
        None means "all llm_exposed tools".
    """
    name: str
    system_prompt: str
    classifier_predicate: Callable[[str], bool]
    max_iterations: int = 8
    critic_step: Optional[CriticStep] = None  # populated in Phase 6
    llm_visible_tools: Optional[List[str]] = None


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

class WorkflowRegistry:
    """
    Ordered registry of WorkflowDefinition objects.

    Registration
    ------------
    Call register(definition) in priority order — first matching predicate wins.
    Use is_fallback=True for the catch-all that handles everything not matched earlier.

    Lookup
    ------
    get_for_task(task_type_str) iterates non-fallback definitions first, then
    falls back to the registered fallback if no predicate matched.  Raises
    RuntimeError if no fallback has been registered.
    """

    def __init__(self) -> None:
        self._definitions: List[WorkflowDefinition] = []
        self._fallback: Optional[WorkflowDefinition] = None

    def register(
        self,
        definition: WorkflowDefinition,
        is_fallback: bool = False,
    ) -> "WorkflowRegistry":
        """Register *definition*. Returns self for chaining."""
        if is_fallback:
            self._fallback = definition
        else:
            self._definitions.append(definition)
        return self

    def get_for_task(self, task_type: str) -> WorkflowDefinition:
        """
        Return the first registered definition whose predicate matches *task_type*.
        Falls back to the fallback definition if no specific one matched.

        Raises RuntimeError if no fallback is registered and no predicate matched.
        """
        for defn in self._definitions:
            if defn.classifier_predicate(task_type):
                return defn
        if self._fallback is not None:
            return self._fallback
        raise RuntimeError(
            f"WorkflowRegistry: no workflow matched task_type={task_type!r} "
            "and no fallback is registered."
        )

    def names(self) -> List[str]:
        """Return all registered workflow names (specific + fallback)."""
        result = [d.name for d in self._definitions]
        if self._fallback is not None:
            result.append(self._fallback.name)
        return result


# ---------------------------------------------------------------------------
# Build the global registry
# ---------------------------------------------------------------------------

def _build_workflow_registry() -> WorkflowRegistry:
    from .retro_prompts import NATIVE_RETRO_SYSTEM_PROMPT
    from .forward_prompts import NATIVE_FORWARD_SYSTEM_PROMPT
    from .prompts import NATIVE_SYSTEM_PROMPT

    _COMPOSITE_SUPPORT_TOOLS = [
        "resolve_chemical",
        "reagent_assistant",
    ]
    _SPECIALIST_ANALYSIS_TOOLS = [
        "featurize_molecule",
        "assess_snar_feasibility",
    ]
    _RETRO_FACADE_AND_SPECIALISTS = [
        "retrosynthesis_step",
        "plan_route",
        "apply_hte_templates",
        "search_by_product_similarity",
        "validate_synthesis_proposal",
    ]
    _FORWARD_FACADE_AND_SPECIALISTS = [
        "forward_synthesis_step",
        "plan_forward_route",
        "validate_synthesis_proposal",
        "featurize_molecule",
    ]
    _GENERAL_FACADE_AND_SPECIALISTS = [
        "analyze_reaction",
        "validate_synthesis_proposal",
        "featurize_molecule",
        "assess_snar_feasibility",
    ]

    registry = WorkflowRegistry()

    registry.register(
        WorkflowDefinition(
            name="retrosynthesis",
            system_prompt=NATIVE_RETRO_SYSTEM_PROMPT,
            classifier_predicate=lambda t: t == "retrosynthesis",
            max_iterations=10,
            critic_step=CriticStep(enabled=True, max_findings=5, min_severity="warning"),
            llm_visible_tools=list(dict.fromkeys(_RETRO_FACADE_AND_SPECIALISTS + _COMPOSITE_SUPPORT_TOOLS)),
        )
    )

    registry.register(
        WorkflowDefinition(
            name="forward_synthesis",
            system_prompt=NATIVE_FORWARD_SYSTEM_PROMPT,
            classifier_predicate=lambda t: t == "forward_synthesis",
            max_iterations=8,
            critic_step=None,
            llm_visible_tools=list(dict.fromkeys(_FORWARD_FACADE_AND_SPECIALISTS + _COMPOSITE_SUPPORT_TOOLS)),
        )
    )

    registry.register(
        WorkflowDefinition(
            name="forward_chemistry",
            system_prompt=NATIVE_SYSTEM_PROMPT,
            classifier_predicate=lambda t: True,   # catch-all fallback
            max_iterations=8,
            critic_step=None,
            llm_visible_tools=list(dict.fromkeys(_GENERAL_FACADE_AND_SPECIALISTS + _COMPOSITE_SUPPORT_TOOLS)),
        ),
        is_fallback=True,
    )

    return registry


# Singleton — imported throughout the codebase
WORKFLOW_REGISTRY: WorkflowRegistry = _build_workflow_registry()
