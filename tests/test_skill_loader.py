from __future__ import annotations

from pathlib import Path

from chem_coworker.skills.loader import SkillLoader
from chem_coworker.skills.manifest import parse_skill_markdown


def _write_skill(root: Path, folder: str, *, skill_id: str, summary: str) -> None:
    target = root / folder
    target.mkdir(parents=True, exist_ok=True)
    (target / "SKILL.md").write_text(
        f"""---
id: {skill_id}
name: Demo Skill
version: 1
summary: {summary}
category: chemistry
workflow_targets:
  - forward_chemistry
tool_allowlist:
  - analyze_reaction
eligibility:
  rdkit: false
  python_modules: []
  data_files: []
  env_vars: []
  taxonomy_ids: []
prompting:
  inject_mode: on_demand
  max_tokens_hint: 10
answer_contract:
  require_tool_evidence: true
priority: 50
---

## Notes

Instruction body.
""",
        encoding="utf-8",
    )


def test_parse_skill_markdown_extracts_manifest_and_body() -> None:
    manifest = parse_skill_markdown(
        """---
id: test_skill
name: Test Skill
version: 2
summary: Sample summary
category: chemistry
workflow_targets:
  - forward_chemistry
tool_allowlist:
  - analyze_reaction
eligibility:
  rdkit: false
  python_modules: []
  data_files: []
  env_vars: []
  taxonomy_ids: []
prompting:
  inject_mode: on_demand
  max_tokens_hint: 12
answer_contract:
  require_tool_evidence: true
priority: 51
---

## Instructions

Use evidence first.
""",
        source_path="demo/SKILL.md",
    )

    assert manifest.id == "test_skill"
    assert manifest.version == 2
    assert manifest.tool_allowlist == ["analyze_reaction"]
    assert "Use evidence first." in manifest.instructions_md


def test_skill_loader_workspace_overrides_user_and_bundled(tmp_path) -> None:
    bundled = tmp_path / "bundled"
    user = tmp_path / "user"
    workspace = tmp_path / "workspace"

    _write_skill(bundled, "same", skill_id="shared", summary="bundled summary")
    _write_skill(user, "same", skill_id="shared", summary="user summary")
    _write_skill(workspace, "same", skill_id="shared", summary="workspace summary")

    loader = SkillLoader(
        bundled_dir=bundled,
        user_dir=user,
        workspace_dir=workspace,
        workspace_root=tmp_path,
    )
    loaded = loader.load_all()

    assert len(loaded) == 1
    assert loaded[0].manifest.id == "shared"
    assert loaded[0].manifest.summary == "workspace summary"
    assert loaded[0].source_label == "workspace"
