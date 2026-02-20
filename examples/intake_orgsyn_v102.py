"""
Batch intake script — Organic Syntheses Volume 102

Reads URLs from orgsyn_v102_valid_p0001_p0600.csv and extracts chemistry notes
for each procedure using the NotesExtractor pipeline.

Usage:
    python examples/intake_orgsyn_v102.py
    python examples/intake_orgsyn_v102.py --model deepseek-v3.2 --provider aliyun
    python examples/intake_orgsyn_v102.py --dry-run          # show URLs without fetching
    python examples/intake_orgsyn_v102.py --resume           # skip already-fetched URLs
    python examples/intake_orgsyn_v102.py --delay 5          # seconds between requests

Results are appended to: examples/intake_orgsyn_v102_log.jsonl
"""
from __future__ import annotations

import argparse
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
# Helpers
# ---------------------------------------------------------------------------

_CSV = Path(__file__).parent / "orgsyn_v102_valid_p0001_p0600.csv"
_LOG = Path(__file__).parent / "intake_orgsyn_v102_log.jsonl"

# ANSI colors (minimal — no external dep)
_R    = "\033[0m"
_BOLD = "\033[1m"
_DIM  = "\033[2m"
_GRN  = "\033[92m"
_YEL  = "\033[93m"
_RED  = "\033[91m"
_CYN  = "\033[96m"


def _read_urls() -> list[str]:
    """Load and clean the URL list from the CSV (one URL per line)."""
    text = _CSV.read_text(encoding="utf-8")
    urls = [line.strip() for line in text.splitlines() if line.strip()]
    return urls


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
        description="Batch intake of Organic Syntheses v102 procedures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--model",    default="deepseek-v3.2", help="LLM model for extraction")
    parser.add_argument("--provider", default="aliyun",        help="LLM provider (openai or aliyun)")
    parser.add_argument("--note-type", default="reactions",
                        choices=["reactions", "mechanisms", "substrates", "protocols"],
                        help="Notes subfolder to write to (default: reactions)")
    parser.add_argument("--delay",    type=float, default=2.0,
                        help="Seconds to wait between requests (default: 2)")
    parser.add_argument("--resume",   action="store_true",
                        help="Skip URLs that already succeeded in a previous run")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print URLs without fetching anything")
    parser.add_argument("--no-save",  action="store_true",
                        help="Do not save fetched pages to literature/ folder")
    parser.add_argument("--verbose",  action="store_true",
                        help="Show full extractor output")
    args = parser.parse_args()

    urls = _read_urls()
    total = len(urls)

    print(f"\n  {_BOLD}◆ OrgSyn v102 Batch Intake{_R}")
    print(f"  {'─' * 44}")
    print(f"  {_DIM}Model  {_R}  {_BOLD}{args.model}{_R}  {_CYN}{args.provider}{_R}")
    print(f"  {_DIM}URLs   {_R}  {total}")
    print(f"  {_DIM}Log    {_R}  {_DIM}{_LOG.name}{_R}")
    print(f"  {_DIM}Delay  {_R}  {args.delay}s between requests")
    if args.resume:
        print(f"  {_YEL}Resume mode ON — skipping previously successful URLs{_R}")
    if args.dry_run:
        print(f"  {_YEL}Dry-run mode — no requests will be made{_R}")
    print()

    if args.dry_run:
        for i, url in enumerate(urls, 1):
            print(f"  {_DIM}{i:>2}.{_R}  {url}")
        print()
        return

    # Check which URLs are already done (resume mode)
    done_urls: set[str] = _already_done(_LOG) if args.resume else set()
    if done_urls:
        print(f"  {_DIM}Skipping {len(done_urls)} already-completed URLs{_R}\n")

    # Initialise extractor once (reuse HTTP session and LLM client)
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

        prefix = f"  {_DIM}{idx:>2}/{total}{_R}"

        if url in done_urls:
            print(f"{prefix}  {_DIM}skip  {url}{_R}")
            n_skip += 1
            continue

        # Progress bar
        bar = _bar(idx - 1, total)
        print(f"{prefix}  {_DIM}{bar}{_R}  {_DIM}{url}{_R}")

        t0 = time.time()
        try:
            result = extractor.intake(
                source=url,
                reaction_type="",          # let LLM auto-detect
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
                _append_log(_LOG, {
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
                _append_log(_LOG, {
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
            _append_log(_LOG, {
                "url": url,
                "success": False,
                "error": str(exc),
                "elapsed_s": round(elapsed, 2),
                "timestamp": datetime.now().isoformat(),
            })
            n_fail += 1

        # Polite delay between requests
        if idx < total and args.delay > 0:
            time.sleep(args.delay)

    # Summary
    print()
    print(f"  {'─' * 44}")
    print(
        f"  {_GRN}✓ {n_ok} succeeded{_R}  "
        + (f"{_DIM}{n_skip} skipped{_R}  " if n_skip else "")
        + (f"{_RED}✗ {n_fail} failed{_R}" if n_fail else "")
    )
    if n_fail:
        print(f"  {_DIM}Re-run with --resume to retry failed URLs{_R}")
    print(f"  {_DIM}Log: {_LOG}{_R}\n")


if __name__ == "__main__":
    main()
