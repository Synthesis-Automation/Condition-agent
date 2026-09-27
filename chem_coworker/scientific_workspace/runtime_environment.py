"""Discover local search tools and prepare an isolated agent process environment."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Mapping


def _ripgrep_candidates(
    environment: Mapping[str, str], home: Path, runtime_executable: str | None,
) -> list[tuple[Path, str]]:
    path_value = next((value for key, value in environment.items() if key.upper() == "PATH"), "")
    found = shutil.which("rg", path=path_value)
    candidates = [(Path(found), "PATH")] if found else []
    if runtime_executable:
        runtime = Path(runtime_executable)
        if runtime.is_file():
            candidates.append((runtime.parent / ("rg.exe" if os.name == "nt" else "rg"), "runtime_bundle"))
    for editor in (".vscode", ".vscode-insiders"):
        extensions = home / editor / "extensions"
        bundled = list(extensions.glob("openai.chatgpt-*/bin/*/rg.exe"))
        bundled.extend(extensions.glob("openai.chatgpt-*/bin/*/rg"))
        candidates.extend((path, "editor_extension") for path in sorted(
            bundled, key=lambda path: path.stat().st_mtime, reverse=True,
        ))
    roots = []
    if environment.get("LOCALAPPDATA"):
        roots.append(Path(environment["LOCALAPPDATA"]) / "Programs")
    if environment.get("ProgramFiles"):
        roots.append(Path(environment["ProgramFiles"]))
    for root in roots:
        for editor in ("Microsoft VS Code", "Microsoft VS Code Insiders"):
            for modules in ("node_modules", "node_modules.asar.unpacked"):
                candidates.append((root / editor / "resources" / "app" / modules
                                   / "@vscode" / "ripgrep" / "bin" / "rg.exe", "editor_bundle"))
    return candidates


def discover_ripgrep(
    *, environment: Mapping[str, str] | None = None, home: Path | None = None,
    runtime_executable: str | None = None,
) -> dict[str, Any]:
    """Verify PATH or installed-editor ripgrep without installing or changing PATH."""
    environment = os.environ if environment is None else environment
    inspected: set[Path] = set()
    for candidate, source in _ripgrep_candidates(environment, home or Path.home(), runtime_executable):
        candidate = candidate.resolve()
        if candidate in inspected or not candidate.is_file():
            continue
        inspected.add(candidate)
        if candidate.suffix.lower() in {".cmd", ".bat", ".ps1"}:
            continue
        try:
            result = subprocess.run(
                [str(candidate), "--version"], capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=5,
                **({"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {}),
            )
        except (OSError, subprocess.SubprocessError):
            continue
        lines = result.stdout.splitlines()
        if result.returncode == 0 and lines and lines[0].startswith("ripgrep "):
            return {
                "status": "available", "path": str(candidate), "version": lines[0],
                "source": source, "child_path_configured": True,
                "sandbox_execution": "not_checked",
            }
    return {
        "status": "unavailable", "path": None, "child_path_configured": False,
        "fallback": "Use Python pathlib/re or PowerShell Get-ChildItem/Select-String.",
    }


def runtime_environment(
    repository: Path, *, search_tool: Mapping[str, Any],
    environment: Mapping[str, str] | None = None,
) -> dict[str, str]:
    """Configure only the child process; preserve the caller's environment and sandbox."""
    child = dict(os.environ if environment is None else environment)
    child["PYTHONPATH"] = str(repository) + os.pathsep + child.get("PYTHONPATH", "")
    child["PYTHONIOENCODING"] = "utf-8"
    child["PYTHONDONTWRITEBYTECODE"] = "1"
    if search_tool.get("status") == "available" and search_tool.get("path"):
        directory = str(Path(search_tool["path"]).parent)
        # Windows environment keys are case-insensitive, but a copied plain dict
        # need not be. Reuse its existing spelling instead of creating two keys.
        path_key = next((key for key in child if key.upper() == "PATH"), "PATH")
        paths = child.get(path_key, "").split(os.pathsep)
        paths = [path for path in paths if os.path.normcase(path) != os.path.normcase(directory)]
        child[path_key] = os.pathsep.join([directory, *paths])
    return child
