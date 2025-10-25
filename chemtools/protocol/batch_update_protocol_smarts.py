"""
Batch updater for protocol SMARTS applicability patterns.

The historical tooling exposed a command-line utility that iterated over the
protocol database, generated SMARTS applicability patterns, and embedded them
back into each JSON file.  This module re-implements the small surface area
required by the test-suite, focusing on deterministic, easy-to-reason-about
behaviour.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

from .smarts_generator_cli import ReactionSmartsApplicability, SmartsGenerator


@dataclass
class ProcessingResult:
    """Outcome of processing a single protocol file."""

    success: bool
    filename: str
    message: str
    old_pattern: Optional[dict] = None


class ProtocolSmartsUpdater:
    """Update ``reaction_smarts_applicability`` entries within protocol JSON files."""

    def __init__(self, protocol_dir: Path | str, dry_run: bool = False):
        self.protocol_dir = Path(protocol_dir)
        self.dry_run = dry_run

    # -- helpers -----------------------------------------------------------------

    def _load_protocol(self, path: Path) -> dict:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    def _save_protocol(self, path: Path, payload: dict) -> None:
        if self.dry_run:
            return
        with path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")

    # -- processing --------------------------------------------------------------

    def process_protocol_file(self, path: Path | str) -> ProcessingResult:
        path = Path(path)
        filename = path.name

        try:
            payload = self._load_protocol(path)
        except json.JSONDecodeError as exc:
            return ProcessingResult(
                success=False,
                filename=filename,
                message=f"JSON parsing error: {exc}",
            )

        reaction = payload.get("reaction", {})
        reaction_smiles = reaction.get("reaction_smiles")
        if not isinstance(reaction_smiles, str) or not reaction_smiles.strip():
            return ProcessingResult(
                success=False,
                filename=filename,
                message="No reaction_smiles found in protocol.",
            )

        try:
            generator = SmartsGenerator(reaction_smiles)
            applicability = generator.build_applicability()
        except Exception as exc:  # pragma: no cover - defensive
            return ProcessingResult(
                success=False,
                filename=filename,
                message=f"SMARTS generation failed: {exc}",
            )

        old_pattern = reaction.get("reaction_smarts_applicability")
        reaction["reaction_smarts_applicability"] = applicability.to_dict()
        payload["reaction"] = reaction

        self._save_protocol(path, payload)

        msg = "Pattern generated and applied"
        if self.dry_run:
            msg += " (dry run)"

        return ProcessingResult(
            success=True,
            filename=filename,
            message=msg,
            old_pattern=old_pattern,
        )

    def process_all_protocols(self) -> List[ProcessingResult]:
        results: List[ProcessingResult] = []
        for path in sorted(self._iter_protocol_files()):
            results.append(self.process_protocol_file(path))
        return results

    def _iter_protocol_files(self) -> Iterable[Path]:
        if not self.protocol_dir.exists():
            return []
        return self.protocol_dir.glob("*.json")


__all__ = ["ProcessingResult", "ProtocolSmartsUpdater"]
