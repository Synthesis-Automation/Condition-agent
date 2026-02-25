"""Output renderers for non-interactive CLI commands."""

from __future__ import annotations

import json
import sys
from typing import Any, Dict, Literal, TextIO

from .style import C
from .ui import TerminalUI

AskOutputFormat = Literal["plain", "json"]
BatchOutputFormat = Literal["plain", "jsonl"]


def render_ask_response(response, output_format: AskOutputFormat, verbose: bool, ui: TerminalUI | None) -> None:
    """Render a single `ask` result."""
    if output_format == "json":
        print(response.to_json())
        return

    assert ui is not None
    ui.render_response(response, verbose=verbose)


def emit_jsonl_record(obj: Dict[str, Any], out_fh: TextIO) -> None:
    """Write one JSONL record."""
    out_fh.write(json.dumps(obj, ensure_ascii=False, default=str) + "\n")


def render_batch_success(
    *,
    idx: int,
    total: int,
    query: str,
    response,
    output_format: BatchOutputFormat,
    verbose: bool,
    ui: TerminalUI | None,
    out_fh: TextIO,
) -> None:
    """Render one successful batch result."""
    if output_format == "jsonl":
        record = response.to_dict()
        record["index"] = idx
        emit_jsonl_record(record, out_fh)
        return

    assert ui is not None
    print(f"\n  {C.LABEL}◆ Batch Item {idx}/{total}{C.R}  {C.DIM}{query[:120]}{C.R}")
    ui.render_response(response, verbose=verbose)


def render_batch_failure(
    *,
    idx: int,
    query: str,
    error: Exception,
    output_format: BatchOutputFormat,
    out_fh: TextIO,
) -> None:
    """Render one failed batch result."""
    if output_format == "jsonl":
        emit_jsonl_record(
            {"index": idx, "query": query, "success": False, "error": str(error)},
            out_fh,
        )
        return

    print(f"\n  {C.ERR}✗{C.R}  Batch item {idx} failed: {error}")


def render_batch_summary(
    *,
    total: int,
    succeeded: int,
    failed: int,
    elapsed_s: float,
    output_format: BatchOutputFormat,
    stream: TextIO | None = None,
) -> None:
    """
    Render batch summary stats.

    For JSONL output, summary is written to stderr by default to avoid corrupting data.
    """
    target = stream or (sys.stderr if output_format == "jsonl" else sys.stdout)
    if output_format == "jsonl":
        target.write(
            json.dumps(
                {
                    "type": "batch_summary",
                    "total": total,
                    "succeeded": succeeded,
                    "failed": failed,
                    "elapsed_s": round(elapsed_s, 2),
                }
            )
            + "\n"
        )
        return

    target.write(
        f"\n  {C.LABEL}◆ Batch Summary{C.R}  "
        f"{C.OK if failed == 0 else C.WARN}{succeeded} ok{C.R} / "
        f"{C.ERR if failed else C.DIM}{failed} failed{C.R} / "
        f"{C.META}{total} total · {elapsed_s:.1f}s{C.R}\n"
    )
