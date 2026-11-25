"""Identify released route inputs and prepare both of their intermediate streams."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

if TYPE_CHECKING:
    from .artifacts import PreprocessingProgress


@dataclass(frozen=True)
class ReleasedSource:
    """One publisher release, its bundle directory and any companion step export."""

    archive: Path
    adapter_id: str
    bundle_relative_path: Path
    step_export: Path | None = None

    def prepare(
        self,
        output_root: Path,
        *,
        force: bool,
        progress_callback: Callable[[PreprocessingProgress], None] | None = None,
        cancel_check: Callable[[], bool] | None = None,
    ) -> dict[str, Any]:
        """Write the route and observation files through the source-specific reader."""
        output = output_root / self.bundle_relative_path
        if self.step_export is not None:
            from .route_release import prepare_route_release

            return prepare_route_release(
                self.archive,
                self.step_export,
                output,
                force=force,
                progress_callback=progress_callback,
                cancel_check=cancel_check,
            )
        from .cbyd import prepare_cbyd_release

        return prepare_cbyd_release(
            self.archive.parent,
            output,
            force=force,
            progress_callback=progress_callback,
            cancel_check=cancel_check,
        )


def released_source(path: Path, raw_root: Path | None = None) -> ReleasedSource | None:
    """Recognize a released source by its archive, once per publisher snapshot."""
    if path.name == "uspto.higher-level.routes.jsonl.gz":
        parent = path.resolve().parent.parent
        prefix = raw_root.resolve() / "routes" if raw_root is not None else None
        relative = (
            parent.relative_to(prefix)
            if prefix and parent.is_relative_to(prefix)
            else Path(parent.name)
        )
        return ReleasedSource(
            path,
            "higher_level_route_source.v1",
            Path("routes") / relative,
            parent / "reactions" / "uspto_original.csv",
        )
    match = re.fullmatch(r"cbyd-(v[0-9]+(?:\.[0-9]+)*)\.tar\.gz", path.name)
    if match and path.parent.name == match[1]:
        return ReleasedSource(
            path,
            "cbyd_release.v1",
            Path("routes/cbyd") / match[1],
        )
    return None
