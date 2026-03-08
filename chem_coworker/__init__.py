"""
ChemCoworker — General-purpose chemistry AI agent.

"Claude Code for chemistry": plan → execute tools in parallel → synthesize.

Quick start:
    from chem_coworker import ChemCoworker

    coworker = ChemCoworker(provider="openai", model="o4-mini")
    response = coworker.run("Recommend conditions for Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
    print(response.answer)

Multi-turn:
    history = []
    resp, history = coworker.chat("What reaction is this?", history)
    resp, history = coworker.chat("Now suggest conditions", history)

Adding a new tool:
    1. Create chem_coworker/tools/my_tools.py
    2. Define ToolPlugin objects, call REGISTRY.register(plugin)
    3. Import the module in tools/__init__.py
    → No changes to agent.py needed
"""
from .agent import ChemCoworker, create_coworker
from .tools import REGISTRY, COWORKER_TOOLS
from .response import ChemResponse
from .classifier import TaskClassifier, TaskType
from .event_bus import EventBus, ChemEvent
from .workflow import WorkflowDefinition, WorkflowRegistry, WORKFLOW_REGISTRY, CriticStep
from .critic import CriticAgent, Finding, Severity, CRITIC_SYSTEM_PROMPT
from .skills import SkillManifest, SkillRegistry, build_default_skill_registry

__all__ = [
    "ChemCoworker",
    "create_coworker",
    "ChemResponse",
    "TaskClassifier",
    "TaskType",
    "REGISTRY",
    "COWORKER_TOOLS",
    "EventBus",
    "ChemEvent",
    "WorkflowDefinition",
    "WorkflowRegistry",
    "WORKFLOW_REGISTRY",
    "CriticStep",
    "CriticAgent",
    "Finding",
    "Severity",
    "CRITIC_SYSTEM_PROMPT",
    "SkillManifest",
    "SkillRegistry",
    "build_default_skill_registry",
]

__version__ = "0.1.0"
