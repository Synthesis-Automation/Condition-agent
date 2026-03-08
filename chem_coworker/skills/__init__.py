"""
ChemCoworker skill infrastructure.
"""

from .catalog import (
    format_skill_catalog,
    format_skill_detail,
    format_skill_instruction_block,
    format_skill_prompt_catalog,
    skill_registry_payload,
)
from .eligibility import SkillEligibilityChecker, SkillEligibilityResult
from .loader import LoadedSkill, SkillLoader
from .manifest import (
    SkillAnswerContract,
    SkillEligibility,
    SkillManifest,
    SkillPrompting,
    SkillTriggers,
    parse_skill_markdown,
)
from .registry import SkillRecord, SkillRegistry, build_default_skill_registry

__all__ = [
    "LoadedSkill",
    "SkillAnswerContract",
    "SkillEligibility",
    "SkillEligibilityChecker",
    "SkillEligibilityResult",
    "SkillLoader",
    "SkillManifest",
    "SkillPrompting",
    "SkillRecord",
    "SkillRegistry",
    "SkillTriggers",
    "build_default_skill_registry",
    "format_skill_catalog",
    "format_skill_detail",
    "format_skill_instruction_block",
    "format_skill_prompt_catalog",
    "parse_skill_markdown",
    "skill_registry_payload",
]
