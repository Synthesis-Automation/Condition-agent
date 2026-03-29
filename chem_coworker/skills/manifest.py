"""
Skill manifest parsing for ChemCoworker.

Skills are markdown files with YAML frontmatter plus free-form instructions.
The frontmatter is structured metadata; the markdown body is only loaded for
inspection or later prompt hydration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class SkillTriggers:
    keywords: list[str] = field(default_factory=list)
    requires_any: list[str] = field(default_factory=list)


@dataclass
class SkillEligibility:
    rdkit: bool = False
    python_modules: list[str] = field(default_factory=list)
    data_files: list[str] = field(default_factory=list)
    env_vars: list[str] = field(default_factory=list)
    taxonomy_ids: list[str] = field(default_factory=list)
    reaction_families: list[str] = field(default_factory=list)
    provider_models: list[str] = field(default_factory=list)


@dataclass
class SkillPrompting:
    inject_mode: str = "on_demand"
    max_tokens_hint: int = 0


@dataclass
class SkillAnswerContract:
    require_tool_evidence: bool = False
    require_taxonomy_alignment: bool = False
    must_surface_warnings: bool = False
    forbid_knowledge_only_conditions: bool = False


@dataclass
class SkillManifest:
    id: str
    name: str
    version: int
    summary: str
    category: str
    tool_policy: str | None = None
    workflow_targets: list[str] = field(default_factory=list)
    triggers: SkillTriggers = field(default_factory=SkillTriggers)
    tool_allowlist: list[str] = field(default_factory=list)
    tool_preferred_order: list[str] = field(default_factory=list)
    tool_default_args: dict[str, dict[str, Any]] = field(default_factory=dict)
    provides_context: list[str] = field(default_factory=list)
    requires_context: list[str] = field(default_factory=list)
    eligibility: SkillEligibility = field(default_factory=SkillEligibility)
    prompting: SkillPrompting = field(default_factory=SkillPrompting)
    answer_contract: SkillAnswerContract = field(default_factory=SkillAnswerContract)
    priority: int = 50
    instructions_md: str = ""
    source_path: str = ""

    @property
    def source_dir(self) -> Path:
        return Path(self.source_path).resolve().parent


def _split_frontmatter(text: str, source_path: str) -> tuple[dict[str, Any], str]:
    if not text.startswith("---"):
        raise ValueError(f"Skill file '{source_path}' is missing YAML frontmatter.")
    lines = text.splitlines()
    end_idx = None
    for idx in range(1, len(lines)):
        if lines[idx].strip() == "---":
            end_idx = idx
            break
    if end_idx is None:
        raise ValueError(f"Skill file '{source_path}' has unterminated YAML frontmatter.")
    frontmatter_text = "\n".join(lines[1:end_idx])
    body = "\n".join(lines[end_idx + 1 :]).strip()
    data = yaml.safe_load(frontmatter_text) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Skill file '{source_path}' frontmatter must decode to a mapping.")
    return data, body


def _as_str(value: Any, field_name: str, source_path: str, *, required: bool = False) -> str:
    if value is None:
        if required:
            raise ValueError(f"Skill file '{source_path}' is missing required field '{field_name}'.")
        return ""
    if not isinstance(value, str):
        raise ValueError(f"Skill file '{source_path}' field '{field_name}' must be a string.")
    return value.strip()


def _as_int(value: Any, field_name: str, source_path: str, *, default: int = 0) -> int:
    if value is None:
        return default
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Skill file '{source_path}' field '{field_name}' must be an integer.")
    return value


def _as_bool(value: Any, field_name: str, source_path: str, *, default: bool = False) -> bool:
    if value is None:
        return default
    if not isinstance(value, bool):
        raise ValueError(f"Skill file '{source_path}' field '{field_name}' must be a boolean.")
    return value


def _as_str_list(value: Any, field_name: str, source_path: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValueError(f"Skill file '{source_path}' field '{field_name}' must be a list.")
    out: list[str] = []
    for item in value:
        if not isinstance(item, str):
            raise ValueError(
                f"Skill file '{source_path}' field '{field_name}' must contain only strings."
            )
        cleaned = item.strip()
        if cleaned:
            out.append(cleaned)
    return out


def parse_skill_markdown(text: str, *, source_path: str = "") -> SkillManifest:
    """Parse one ChemCoworker skill manifest markdown file into a typed manifest."""
    data, instructions_md = _split_frontmatter(text, source_path)

    triggers_data = data.get("triggers") or {}
    if not isinstance(triggers_data, dict):
        raise ValueError(f"Skill file '{source_path}' field 'triggers' must be a mapping.")

    eligibility_data = data.get("eligibility") or {}
    if not isinstance(eligibility_data, dict):
        raise ValueError(f"Skill file '{source_path}' field 'eligibility' must be a mapping.")

    prompting_data = data.get("prompting") or {}
    if not isinstance(prompting_data, dict):
        raise ValueError(f"Skill file '{source_path}' field 'prompting' must be a mapping.")

    tool_default_args_raw = data.get("tool_default_args") or {}
    if not isinstance(tool_default_args_raw, dict):
        raise ValueError(f"Skill file '{source_path}' field 'tool_default_args' must be a mapping.")
    tool_default_args: dict[str, dict[str, Any]] = {}
    for tool_name, args in tool_default_args_raw.items():
        if not isinstance(args, dict):
            raise ValueError(
                f"Skill file '{source_path}' field 'tool_default_args.{tool_name}' must be a mapping."
            )
        tool_default_args[str(tool_name)] = dict(args)

    contract_data = data.get("answer_contract") or {}
    if not isinstance(contract_data, dict):
        raise ValueError(f"Skill file '{source_path}' field 'answer_contract' must be a mapping.")

    manifest = SkillManifest(
        id=_as_str(data.get("id"), "id", source_path, required=True),
        name=_as_str(data.get("name"), "name", source_path, required=True),
        version=_as_int(data.get("version"), "version", source_path, default=1),
        summary=_as_str(data.get("summary"), "summary", source_path, required=True),
        category=_as_str(data.get("category"), "category", source_path, required=True),
        tool_policy=_as_str(data.get("tool_policy"), "tool_policy", source_path) or None,
        workflow_targets=_as_str_list(data.get("workflow_targets"), "workflow_targets", source_path),
        triggers=SkillTriggers(
            keywords=_as_str_list(triggers_data.get("keywords"), "triggers.keywords", source_path),
            requires_any=_as_str_list(
                triggers_data.get("requires_any"), "triggers.requires_any", source_path
            ),
        ),
        tool_allowlist=_as_str_list(data.get("tool_allowlist"), "tool_allowlist", source_path),
        tool_preferred_order=_as_str_list(
            data.get("tool_preferred_order"), "tool_preferred_order", source_path
        ),
        tool_default_args=tool_default_args,
        provides_context=_as_str_list(data.get("provides_context"), "provides_context", source_path),
        requires_context=_as_str_list(data.get("requires_context"), "requires_context", source_path),
        eligibility=SkillEligibility(
            rdkit=_as_bool(eligibility_data.get("rdkit"), "eligibility.rdkit", source_path),
            python_modules=_as_str_list(
                eligibility_data.get("python_modules"), "eligibility.python_modules", source_path
            ),
            data_files=_as_str_list(
                eligibility_data.get("data_files"), "eligibility.data_files", source_path
            ),
            env_vars=_as_str_list(
                eligibility_data.get("env_vars"), "eligibility.env_vars", source_path
            ),
            taxonomy_ids=_as_str_list(
                eligibility_data.get("taxonomy_ids"), "eligibility.taxonomy_ids", source_path
            ),
            reaction_families=_as_str_list(
                eligibility_data.get("reaction_families"),
                "eligibility.reaction_families",
                source_path,
            ),
            provider_models=_as_str_list(
                eligibility_data.get("provider_models"),
                "eligibility.provider_models",
                source_path,
            ),
        ),
        prompting=SkillPrompting(
            inject_mode=_as_str(prompting_data.get("inject_mode"), "prompting.inject_mode", source_path)
            or "on_demand",
            max_tokens_hint=_as_int(
                prompting_data.get("max_tokens_hint"),
                "prompting.max_tokens_hint",
                source_path,
                default=0,
            ),
        ),
        answer_contract=SkillAnswerContract(
            require_tool_evidence=_as_bool(
                contract_data.get("require_tool_evidence"),
                "answer_contract.require_tool_evidence",
                source_path,
            ),
            require_taxonomy_alignment=_as_bool(
                contract_data.get("require_taxonomy_alignment"),
                "answer_contract.require_taxonomy_alignment",
                source_path,
            ),
            must_surface_warnings=_as_bool(
                contract_data.get("must_surface_warnings"),
                "answer_contract.must_surface_warnings",
                source_path,
            ),
            forbid_knowledge_only_conditions=_as_bool(
                contract_data.get("forbid_knowledge_only_conditions"),
                "answer_contract.forbid_knowledge_only_conditions",
                source_path,
            ),
        ),
        priority=_as_int(data.get("priority"), "priority", source_path, default=50),
        instructions_md=instructions_md,
        source_path=str(Path(source_path)) if source_path else "",
    )
    if not manifest.id:
        raise ValueError(f"Skill file '{source_path}' field 'id' cannot be empty.")
    if not manifest.name:
        raise ValueError(f"Skill file '{source_path}' field 'name' cannot be empty.")
    return manifest

