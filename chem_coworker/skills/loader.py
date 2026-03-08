"""
Filesystem loader for ChemCoworker skills.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .manifest import SkillManifest, parse_skill_markdown


@dataclass
class LoadedSkill:
    manifest: SkillManifest
    source_label: str


class SkillLoader:
    def __init__(
        self,
        *,
        bundled_dir: Path | None = None,
        user_dir: Path | None = None,
        workspace_dir: Path | None = None,
        workspace_root: Path | None = None,
    ) -> None:
        repo_root = Path(__file__).resolve().parents[2]
        self.workspace_root = Path(workspace_root).resolve() if workspace_root else Path.cwd().resolve()
        self.bundled_dir = bundled_dir or (repo_root / "chem_coworker_skillpacks" / "bundled")
        self.user_dir = user_dir or (Path.home() / ".chem_coworker" / "skills")
        self.workspace_dir = workspace_dir or (self.workspace_root / "skills")

    def source_directories(self) -> list[tuple[str, Path]]:
        return [
            ("bundled", Path(self.bundled_dir)),
            ("user", Path(self.user_dir)),
            ("workspace", Path(self.workspace_dir)),
        ]

    def _iter_skill_files(self, root: Path) -> list[Path]:
        if not root.exists() or not root.is_dir():
            return []
        files: list[Path] = []
        for child in sorted(root.iterdir(), key=lambda p: p.name.lower()):
            if child.is_dir():
                skill_file = child / "SKILL.md"
                if skill_file.exists():
                    files.append(skill_file)
        return files

    def load_all(self) -> list[LoadedSkill]:
        by_id: dict[str, LoadedSkill] = {}
        for source_label, root in self.source_directories():
            for skill_file in self._iter_skill_files(root):
                text = skill_file.read_text(encoding="utf-8")
                manifest = parse_skill_markdown(text, source_path=str(skill_file))
                by_id[manifest.id] = LoadedSkill(manifest=manifest, source_label=source_label)
        return list(by_id.values())

