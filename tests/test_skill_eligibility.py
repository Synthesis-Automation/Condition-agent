from __future__ import annotations

from pathlib import Path

from chem_coworker.skills.eligibility import SkillEligibilityChecker
from chem_coworker.skills.manifest import SkillEligibility, SkillManifest


def _make_manifest(tmp_path: Path) -> SkillManifest:
    return SkillManifest(
        id="eligibility_test",
        name="Eligibility Test",
        version=1,
        summary="demo",
        category="chemistry",
        eligibility=SkillEligibility(
            rdkit=False,
            python_modules=[],
            data_files=[],
            env_vars=[],
            taxonomy_ids=[],
        ),
        source_path=str(tmp_path / "skills" / "eligibility_test" / "SKILL.md"),
    )


def test_eligibility_checker_reports_missing_env_and_data_file(monkeypatch, tmp_path) -> None:
    manifest = _make_manifest(tmp_path)
    manifest.eligibility.env_vars = ["MISSING_SKILL_ENV"]
    manifest.eligibility.data_files = ["data/not_here"]
    monkeypatch.delenv("MISSING_SKILL_ENV", raising=False)

    checker = SkillEligibilityChecker(workspace_root=tmp_path)
    result = checker.check(manifest)

    assert result.eligible is False
    assert "missing environment variable: MISSING_SKILL_ENV" in result.reasons
    assert "missing data file: data/not_here" in result.reasons


def test_eligibility_checker_accepts_present_env_and_data_file(monkeypatch, tmp_path) -> None:
    manifest = _make_manifest(tmp_path)
    manifest.eligibility.env_vars = ["PRESENT_SKILL_ENV"]
    manifest.eligibility.data_files = ["data/demo"]
    monkeypatch.setenv("PRESENT_SKILL_ENV", "ok")
    data_dir = tmp_path / "data" / "demo"
    data_dir.mkdir(parents=True)

    checker = SkillEligibilityChecker(workspace_root=tmp_path)
    result = checker.check(manifest)

    assert result.eligible is True
    assert result.reasons == []
