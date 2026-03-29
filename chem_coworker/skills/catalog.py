"""
Plain formatting helpers for skill listings.
"""

from __future__ import annotations

from .registry import SkillRecord, SkillRegistry


def format_skill_catalog(records: list[SkillRecord]) -> str:
    if not records:
        return "No eligible skills found."
    lines: list[str] = []
    for record in records:
        manifest = record.manifest
        lines.append(f"- {manifest.id}: {manifest.summary}")
        workflow_text = ", ".join(manifest.workflow_targets) if manifest.workflow_targets else "all"
        tools_text = ", ".join(manifest.tool_allowlist) if manifest.tool_allowlist else "(none declared)"
        lines.append(f"  Workflows: {workflow_text}")
        lines.append(f"  Tools: {tools_text}")
    return "\n".join(lines)


def format_skill_prompt_catalog(records: list[SkillRecord]) -> str:
    """Compact metadata-only block suitable for prompt injection."""
    if not records:
        return "Available skills: (none)"
    lines = ["Available skills:"]
    for record in records:
        manifest = record.manifest
        tools_text = ", ".join(manifest.tool_allowlist[:4]) if manifest.tool_allowlist else "(none declared)"
        lines.append(f"- {manifest.id}: {manifest.summary}")
        lines.append(f"  Tools: {tools_text}")
    return "\n".join(lines)


def format_skill_instruction_block(records: list[SkillRecord]) -> str:
    """Full instruction block for the active skill set."""
    if not records:
        return ""
    lines = ["Active skill instructions:"]
    for record in records:
        manifest = record.manifest
        lines.append(f"## {manifest.name} [{manifest.id}]")
        lines.append(manifest.instructions_md.strip() or "(no additional instructions)")
        # Surface answer_contract constraints so the LLM knows the rules
        contract = manifest.answer_contract
        contract_lines: list[str] = []
        if contract.require_tool_evidence:
            contract_lines.append("- All claims MUST be backed by tool output evidence.")
        if contract.require_taxonomy_alignment:
            contract_lines.append("- Confirm taxonomy-backed reaction identity before recommending.")
        if contract.forbid_knowledge_only_conditions:
            contract_lines.append("- Do NOT fabricate conditions from general knowledge; only report tool-returned results.")
        if contract.must_surface_warnings:
            contract_lines.append("- Surface all tool warnings and caveats explicitly in the answer.")
        if contract_lines:
            lines.append("\n**Answer contract (enforced):**")
            lines.extend(contract_lines)
        # Surface tool_default_args so the LLM uses the right parameters
        if manifest.tool_default_args:
            default_lines: list[str] = []
            for tool_name, defaults in manifest.tool_default_args.items():
                parts = ", ".join(f"{k}={v!r}" for k, v in defaults.items())
                default_lines.append(f"- {tool_name}: use {parts}")
            if default_lines:
                lines.append("\n**Recommended tool parameters:**")
                lines.extend(default_lines)
    return "\n\n".join(lines)


def format_skill_detail(record: SkillRecord) -> str:
    manifest = record.manifest
    lines = [
        f"ID: {manifest.id}",
        f"Name: {manifest.name}",
        f"Version: {manifest.version}",
        f"Source: {record.source_label}",
        f"Eligible: {'yes' if record.eligibility.eligible else 'no'}",
        f"Summary: {manifest.summary}",
        f"Category: {manifest.category}",
    ]
    if manifest.workflow_targets:
        lines.append(f"Workflows: {', '.join(manifest.workflow_targets)}")
    if manifest.tool_policy:
        lines.append(f"Tool policy: {manifest.tool_policy}")
    if manifest.tool_allowlist:
        lines.append(f"Tool allowlist: {', '.join(manifest.tool_allowlist)}")
    if record.eligibility.reasons:
        lines.append(f"Eligibility notes: {'; '.join(record.eligibility.reasons)}")
    if manifest.instructions_md:
        lines.append("")
        lines.append(manifest.instructions_md)
    return "\n".join(lines)


def skill_registry_payload(registry: SkillRegistry) -> dict:
    eligible = registry.eligible_records()
    suppressed = registry.suppressed_records()
    return {
        "eligible_count": len(eligible),
        "suppressed_count": len(suppressed),
        "eligible": [
            {
                "id": record.manifest.id,
                "name": record.manifest.name,
                "summary": record.manifest.summary,
                "source": record.source_label,
                "workflow_targets": list(record.manifest.workflow_targets),
                "tool_allowlist": list(record.manifest.tool_allowlist),
                "priority": record.manifest.priority,
            }
            for record in eligible
        ],
        "suppressed": [
            {
                "id": record.manifest.id,
                "name": record.manifest.name,
                "source": record.source_label,
                "reasons": list(record.eligibility.reasons),
            }
            for record in suppressed
        ],
    }
