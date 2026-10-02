"""Installed workspace resources and the repository containing domain packages."""

from pathlib import Path

WORKSPACE_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = WORKSPACE_ROOT.parent.parent
