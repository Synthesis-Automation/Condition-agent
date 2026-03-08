"""
Eligibility checks for skills.
"""

from __future__ import annotations

import importlib.util
import os
from dataclasses import dataclass, field
from pathlib import Path

from .manifest import SkillManifest


@dataclass
class SkillEligibilityResult:
    skill_id: str
    eligible: bool
    reasons: list[str] = field(default_factory=list)


class SkillEligibilityChecker:
    def __init__(self, workspace_root: Path | None = None) -> None:
        self.workspace_root = Path(workspace_root).resolve() if workspace_root else Path.cwd().resolve()
        self.repo_root = Path(__file__).resolve().parents[2]

    def _candidate_paths(self, raw_path: str) -> list[Path]:
        path = Path(raw_path)
        if path.is_absolute():
            return [path]
        return [
            self.workspace_root / path,
            self.repo_root / path,
        ]

    def _check_taxonomy_id(self, taxonomy_id: str) -> bool:
        try:
            from chemtools.taxonomy.reaction_catalog import get_reaction_type
        except Exception:
            return False
        try:
            return get_reaction_type(taxonomy_id) is not None
        except Exception:
            return False

    def check(self, manifest: SkillManifest) -> SkillEligibilityResult:
        reasons: list[str] = []
        eligibility = manifest.eligibility

        if eligibility.rdkit and importlib.util.find_spec("rdkit") is None:
            reasons.append("missing python module: rdkit")

        for module_name in eligibility.python_modules:
            if importlib.util.find_spec(module_name) is None:
                reasons.append(f"missing python module: {module_name}")

        for raw_path in eligibility.data_files:
            candidates = self._candidate_paths(raw_path)
            if not any(candidate.exists() for candidate in candidates):
                reasons.append(f"missing data file: {raw_path}")

        for env_var in eligibility.env_vars:
            if not os.getenv(env_var):
                reasons.append(f"missing environment variable: {env_var}")

        for taxonomy_id in eligibility.taxonomy_ids:
            if not self._check_taxonomy_id(taxonomy_id):
                reasons.append(f"unknown taxonomy id: {taxonomy_id}")

        return SkillEligibilityResult(
            skill_id=manifest.id,
            eligible=not reasons,
            reasons=reasons,
        )

