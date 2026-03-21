"""
Registry for loaded ChemCoworker skills.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .eligibility import SkillEligibilityChecker, SkillEligibilityResult
from .loader import SkillLoader
from .manifest import SkillManifest


@dataclass
class SkillRecord:
    manifest: SkillManifest
    source_label: str
    eligibility: SkillEligibilityResult


class SkillRegistry:
    def __init__(self, records: list[SkillRecord]) -> None:
        self._records = list(records)
        self._by_id = {record.manifest.id: record for record in records}

    def get(self, skill_id: str) -> SkillManifest | None:
        record = self._by_id.get(skill_id)
        return record.manifest if record else None

    def get_record(self, skill_id: str) -> SkillRecord | None:
        return self._by_id.get(skill_id)

    def all_records(self) -> list[SkillRecord]:
        return list(self._records)

    def eligible_records(self) -> list[SkillRecord]:
        return sorted(
            (record for record in self._records if record.eligibility.eligible),
            key=lambda record: (-record.manifest.priority, record.manifest.id),
        )

    def suppressed_records(self) -> list[SkillRecord]:
        return sorted(
            (record for record in self._records if not record.eligibility.eligible),
            key=lambda record: record.manifest.id,
        )

    def catalog_for_workflow(self, workflow_name: str) -> list[SkillManifest]:
        manifests: list[SkillManifest] = []
        for record in self.eligible_records():
            if not record.manifest.workflow_targets or workflow_name in record.manifest.workflow_targets:
                manifests.append(record.manifest)
        return manifests

    def match_for_query(
        self,
        query: str,
        task_type: str,
        smiles_present: bool,
        workflow_name: str | None = None,
    ) -> list[SkillManifest]:
        query_lc = str(query or "").lower()
        matched: list[SkillManifest] = []
        for record in self.eligible_records():
            manifest = record.manifest
            if manifest.workflow_targets and workflow_name not in set(manifest.workflow_targets):
                continue
            requires_any_ok = not manifest.triggers.requires_any or (
                smiles_present and any(
                    item in {"reaction_smiles", "reaction_context", "molecule_smiles"}
                    for item in manifest.triggers.requires_any
                )
            )
            if not requires_any_ok:
                continue
            keyword_match = any(keyword.lower() in query_lc for keyword in manifest.triggers.keywords)
            if manifest.triggers.keywords:
                if keyword_match:
                    matched.append(manifest)
            else:
                matched.append(manifest)
        return matched


def build_default_skill_registry(workspace_root: Path | None = None) -> SkillRegistry:
    loader = SkillLoader(workspace_root=workspace_root)
    checker = SkillEligibilityChecker(workspace_root=workspace_root)
    records: list[SkillRecord] = []
    for loaded in loader.load_all():
        result = checker.check(loaded.manifest)
        records.append(
            SkillRecord(
                manifest=loaded.manifest,
                source_label=loaded.source_label,
                eligibility=result,
            )
        )
    return SkillRegistry(records)
