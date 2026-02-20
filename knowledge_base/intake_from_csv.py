"""
Batch intake script — extract chemistry notes from a list of URLs in a CSV file.

Reads URLs from a CSV file and runs each through the NotesExtractor pipeline,
saving notes to knowledge_base/notes/ and raw sources to knowledge_base/sources/.

CSV formats supported:
  • One URL per line (no header):
        https://orgsyn.org/demo.aspx?prep=v102p0001
        https://orgsyn.org/demo.aspx?prep=v102p0019

  • CSV with a header column named "url", "URL", "link", or "href":
        url,title,notes
        https://example.com/paper1,Paper 1,

The log file is written alongside the CSV: {csv_stem}_log.jsonl

Usage:
    python knowledge_base/intake_from_csv.py knowledge_base/orgsyn_v102.csv
    python knowledge_base/intake_from_csv.py my_urls.csv --model deepseek-v3.2 --provider aliyun
    python knowledge_base/intake_from_csv.py my_urls.csv --dry-run
    python knowledge_base/intake_from_csv.py my_urls.csv --resume
    python knowledge_base/intake_from_csv.py my_urls.csv --delay 5 --url-column link
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

# Make sure project root is on the path
_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

from dotenv import load_dotenv
load_dotenv()


# ---------------------------------------------------------------------------
# ANSI colors
# ---------------------------------------------------------------------------

_R    = "\033[0m"
_BOLD = "\033[1m"
_DIM  = "\033[2m"
_GRN  = "\033[92m"
_YEL  = "\033[93m"
_RED  = "\033[91m"
_CYN  = "\033[96m"

# Column names recognised as URL columns (case-insensitive, checked in order)
_URL_COLUMNS = ["url", "link", "href", "uri", "address"]


# ---------------------------------------------------------------------------
# CSV reading
# ---------------------------------------------------------------------------

def _read_urls(csv_path: Path, url_column: str = "") -> list[str]:
    """
    Read URLs from a CSV file.

    Supports:
      - Plain list: one URL per line, no header
      - CSV with header: auto-detects a URL column, or uses --url-column
    """
    text = csv_path.read_text(encoding="utf-8-sig")  # strip BOM if present
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        return []

    # --- Detect if first line is a header ---
    first = lines[0]
    # If the first line starts with http, treat the whole file as plain URL list
    if first.lower().startswith("http"):
        return [l for l in lines if l.lower().startswith("http")]

    # Try to parse as CSV with header
    reader = csv.DictReader(text.splitlines())
    fieldnames = [f.strip() for f in (reader.fieldnames or [])]

    # Find URL column
    col = None
    if url_column:
        # Exact match first, then case-insensitive
        if url_column in fieldnames:
            col = url_column
        else:
            for f in fieldnames:
                if f.lower() == url_column.lower():
                    col = f
                    break
        if col is None:
            print(f"{_RED}Error: column '{url_column}' not found. Available: {fieldnames}{_R}")
            sys.exit(1)
    else:
        for candidate in _URL_COLUMNS:
            for f in fieldnames:
                if f.lower() == candidate:
                    col = f
                    break
            if col:
                break

    if col:
        urls = []
        for row in csv.DictReader(text.splitlines()):
            val = (row.get(col) or "").strip()
            if val.lower().startswith("http"):
                urls.append(val)
        return urls

    # Fallback: treat every non-empty line that looks like a URL as a URL
    return [l for l in lines if l.lower().startswith("http")]


# ---------------------------------------------------------------------------
# Resume / log helpers
# ---------------------------------------------------------------------------

def _already_done(log_path: Path) -> set[str]:
    """Return URLs that completed successfully in a previous run."""
    done: set[str] = set()
    if not log_path.exists():
        return done
    for line in log_path.read_text(encoding="utf-8").splitlines():
        try:
            rec = json.loads(line)
            if rec.get("success"):
                done.add(rec["url"])
        except Exception:
            pass
    return done


def _append_log(log_path: Path, record: dict) -> None:
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, default=str) + "\n")


def _bar(done: int, total: int, width: int = 30) -> str:
    filled = int(width * done / total) if total else 0
    return f"[{'█' * filled}{'░' * (width - filled)}]"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch intake of chemistry URLs from a CSV file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "csv",
        type=Path,
        help="Path to the CSV file containing URLs (one per line, or with a URL column)",
    )
    parser.add_argument("--model",    default="deepseek-v3.2", help="LLM model for extraction (default: deepseek-v3.2)")
    parser.add_argument("--provider", default="aliyun",        help="LLM provider: openai or aliyun (default: aliyun)")
    parser.add_argument("--note-type", default="reactions",
                        choices=["reactions", "mechanisms", "substrates", "protocols"],
                        help="Notes subfolder to write to (default: reactions)")
    parser.add_argument("--url-column", default="",
                        help="CSV column name that contains URLs (auto-detected if omitted)")
    parser.add_argument("--delay",    type=float, default=2.0,
                        help="Seconds to wait between requests (default: 2)")
    parser.add_argument("--resume",   action="store_true",
                        help="Skip URLs that already succeeded in a previous run")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print URLs without fetching anything")
    parser.add_argument("--no-save",  action="store_true",
                        help="Do not save fetched pages to knowledge_base/sources/")
    parser.add_argument("--log",      type=Path, default=None,
                        help="Path for the JSONL log file (default: {csv_stem}_log.jsonl alongside the CSV)")
    parser.add_argument("--verbose",  action="store_true",
                        help="Show full extractor output")
    args = parser.parse_args()

    csv_path: Path = args.csv.resolve()
    if not csv_path.exists():
        print(f"{_RED}Error: CSV file not found: {csv_path}{_R}")
        sys.exit(1)

    log_path: Path = args.log or csv_path.with_name(csv_path.stem + "_log.jsonl")

    urls = _read_urls(csv_path, args.url_column)
    total = len(urls)

    if not urls:
        print(f"{_YEL}No URLs found in {csv_path.name}{_R}")
        sys.exit(0)

    print(f"\n  {_BOLD}◆ Batch Intake — {csv_path.name}{_R}")
    print(f"  {'─' * 48}")
    print(f"  {_DIM}CSV    {_R}  {csv_path.name}")
    print(f"  {_DIM}Model  {_R}  {_BOLD}{args.model}{_R}  {_CYN}{args.provider}{_R}")
    print(f"  {_DIM}Type   {_R}  {args.note_type}")
    print(f"  {_DIM}URLs   {_R}  {total}")
    print(f"  {_DIM}Log    {_R}  {_DIM}{log_path.name}{_R}")
    print(f"  {_DIM}Delay  {_R}  {args.delay}s between requests")
    if args.resume:
        print(f"  {_YEL}Resume mode ON — skipping previously successful URLs{_R}")
    if args.dry_run:
        print(f"  {_YEL}Dry-run mode — no requests will be made{_R}")
    print()

    if args.dry_run:
        for i, url in enumerate(urls, 1):
            print(f"  {_DIM}{i:>3}.{_R}  {url}")
        print()
        return

    # Resume: load already-done URLs
    done_urls: set[str] = _already_done(log_path) if args.resume else set()
    if done_urls:
        print(f"  {_DIM}Skipping {len(done_urls)} already-completed URLs{_R}\n")

    # Initialise extractor once (reuses HTTP session and LLM client)
    from chem_coworker.extractor import NotesExtractor
    extractor = NotesExtractor(
        provider=args.provider,
        model=args.model,
        verbose=args.verbose,
    )

    n_ok = n_skip = n_fail = 0

    for idx, url in enumerate(urls, 1):
        url = url.strip()
        if not url:
            continue

        prefix = f"  {_DIM}{idx:>3}/{total}{_R}"

        if url in done_urls:
            print(f"{prefix}  {_DIM}skip  {url}{_R}")
            n_skip += 1
            continue

        bar = _bar(idx - 1, total)
        print(f"{prefix}  {_DIM}{bar}{_R}  {_DIM}{url}{_R}")

        t0 = time.time()
        try:
            result = extractor.intake(
                source=url,
                reaction_type="",
                note_type=args.note_type,
                save_to_literature=not args.no_save,
            )
            elapsed = time.time() - t0

            if result.get("success"):
                types_str = ", ".join(result.get("reaction_types", []))
                chars = result.get("char_count", 0)
                print(
                    f"       {_GRN}✓{_R}  "
                    f"{_DIM}{elapsed:.1f}s{_R}  "
                    f"{chars:,} chars  "
                    f"{_CYN}{types_str}{_R}"
                )
                _append_log(log_path, {
                    "url": url,
                    "success": True,
                    "reaction_types": result.get("reaction_types", []),
                    "notes_files": result.get("notes_files", []),
                    "char_count": chars,
                    "elapsed_s": round(elapsed, 2),
                    "timestamp": datetime.now().isoformat(),
                    "model": args.model,
                    "warnings": result.get("warnings", []),
                })
                n_ok += 1
            else:
                err = result.get("error", "unknown error")
                print(f"       {_RED}✗{_R}  {err}")
                _append_log(log_path, {
                    "url": url,
                    "success": False,
                    "error": err,
                    "elapsed_s": round(elapsed, 2),
                    "timestamp": datetime.now().isoformat(),
                })
                n_fail += 1

        except KeyboardInterrupt:
            print(f"\n  {_YEL}Interrupted — {n_ok} succeeded, {n_fail} failed so far.{_R}")
            print(f"  {_DIM}Run with --resume to continue from here.{_R}\n")
            sys.exit(0)
        except Exception as exc:
            elapsed = time.time() - t0
            print(f"       {_RED}✗{_R}  Exception: {exc}")
            _append_log(log_path, {
                "url": url,
                "success": False,
                "error": str(exc),
                "elapsed_s": round(elapsed, 2),
                "timestamp": datetime.now().isoformat(),
            })
            n_fail += 1

        if idx < total and args.delay > 0:
            time.sleep(args.delay)

    # Summary
    print()
    print(f"  {'─' * 48}")
    print(
        f"  {_GRN}✓ {n_ok} succeeded{_R}  "
        + (f"{_DIM}{n_skip} skipped{_R}  " if n_skip else "")
        + (f"{_RED}✗ {n_fail} failed{_R}" if n_fail else "")
    )
    if n_fail:
        print(f"  {_DIM}Re-run with --resume to retry failed URLs{_R}")
    print(f"  {_DIM}Log: {log_path}{_R}\n")


if __name__ == "__main__":
    main()
