"""
ChemCoworker CLI

Usage:
    python chem_coworker/cli.py                                   # interactive chat
    python chem_coworker/cli.py --model gpt-5.2 --provider openai
    python chem_coworker/cli.py --model glm-5 --provider aliyun
    python chem_coworker/cli.py setup                             # choose & save default model
    python chem_coworker/cli.py intake <source>                   # intake a document
    python chem_coworker/cli.py intake https://...                # fetch + extract URL
    python chem_coworker/cli.py intake paper.pdf --reaction-type suzuki_miyaura

Features:
  - Saved model preference (~/.chemcoworker/config.json)
  - Multi-turn conversation with history
  - Shows plan, tools, and answer in Claude Code-style output
  - Spinner animation during LLM calls
  - Wall-clock timing display
  - Intake subcommand: extract chemistry notes from URLs, PDFs, and text files
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional

if TYPE_CHECKING:
    from chem_coworker.response import ChemResponse


# ---------------------------------------------------------------------------
# ANSI color / style constants
# ---------------------------------------------------------------------------

class C:
    """ANSI color codes. All named for their semantic role."""
    R       = "\033[0m"       # reset
    BOLD    = "\033[1m"
    DIM     = "\033[2m"
    ITALIC  = "\033[3m"

    # Foreground colors (bright variants for readability on dark/light)
    RED     = "\033[91m"
    GREEN   = "\033[92m"
    YELLOW  = "\033[93m"
    BLUE    = "\033[94m"
    MAGENTA = "\033[95m"
    CYAN    = "\033[96m"
    WHITE   = "\033[97m"

    # Semantic aliases — mirrors Claude Code palette
    PROMPT   = "\033[1;94m"   # bold blue  — user prompt chevron
    LABEL    = "\033[1;96m"   # bold cyan  — section labels (◆ Hypothesis, ⎿ Tools …)
    META     = "\033[2m"      # dim        — timing, model, token counts
    TOOL     = "\033[93m"     # yellow     — tool names
    ANSWER   = "\033[97m"     # bright white — answer body
    HYPO     = "\033[96m"     # cyan       — hypothesis text
    CONF     = "\033[2;92m"   # dim green  — confidence score
    WARN     = "\033[91m"     # red        — warnings
    OK       = "\033[92m"     # green      — success marks
    ERR      = "\033[91m"     # red        — error marks
    SECTION  = "\033[2m"      # dim        — separators / structural chrome


# Width for separators
_W = 60
_SEP     = f"{C.SECTION}{'─' * _W}{C.R}"
_SEP_FAT = f"{C.SECTION}{'━' * _W}{C.R}"


# ---------------------------------------------------------------------------
# Model selector
# ---------------------------------------------------------------------------

SELECTABLE_MODELS: List[Dict[str, str]] = [
    {"name": "o4-mini",        "provider": "openai"},   # 1 — default
    {"name": "gpt-5.2",        "provider": "openai"},   # 2
    {"name": "glm-5",          "provider": "aliyun"},   # 3
    {"name": "glm-4.7",        "provider": "aliyun"},   # 4
    {"name": "MiniMax-M2.1",   "provider": "aliyun"},   # 5
    {"name": "deepseek-v3.2",  "provider": "aliyun"},   # 6
]

_ALIYUN_MODELS = {m["name"] for m in SELECTABLE_MODELS if m["provider"] == "aliyun"}


def select_model_interactive() -> Dict[str, str]:
    """Show numbered model menu and return selected {name, provider}."""
    print(f"\n  {C.LABEL}◆ Select model{C.R}")
    print(f"  {C.SECTION}{'─' * 44}{C.R}")
    for i, m in enumerate(SELECTABLE_MODELS, 1):
        default_tag = f"  {C.DIM}← default{C.R}" if i == 1 else ""
        provider_color = C.CYAN if m["provider"] == "openai" else C.MAGENTA
        print(
            f"  {C.DIM}{i}.{C.R}  "
            f"{C.BOLD}{m['name']:<20}{C.R}"
            f"{provider_color}{m['provider']}{C.R}"
            f"{default_tag}"
        )
    print(f"  {C.SECTION}{'─' * 44}{C.R}")
    raw = input(f"  {C.PROMPT}>{C.R} Choice (1–{len(SELECTABLE_MODELS)}, Enter=default): ").strip()
    if not raw:
        return SELECTABLE_MODELS[0]
    try:
        idx = int(raw) - 1
        if 0 <= idx < len(SELECTABLE_MODELS):
            return SELECTABLE_MODELS[idx]
    except ValueError:
        pass
    print(f"  {C.DIM}Invalid — using default (o4-mini){C.R}")
    return SELECTABLE_MODELS[0]


# ---------------------------------------------------------------------------
# Config file (saved model preference)
# ---------------------------------------------------------------------------

_CONFIG_PATH = Path.home() / ".chemcoworker" / "config.json"
_DEFAULT_MODEL = {"name": "o4-mini", "provider": "openai"}


def _load_config() -> Dict[str, str]:
    """Load saved model config. Returns default if file is missing or corrupt."""
    try:
        data = json.loads(_CONFIG_PATH.read_text(encoding="utf-8"))
        if "name" in data and "provider" in data:
            return {"name": data["name"], "provider": data["provider"]}
    except Exception:
        pass
    return _DEFAULT_MODEL.copy()


def _save_config(model: str, provider: str) -> None:
    """Persist model choice to ~/.chemcoworker/config.json."""
    _CONFIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    _CONFIG_PATH.write_text(
        json.dumps({"name": model, "provider": provider}, indent=2),
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Spinner
# ---------------------------------------------------------------------------

class Spinner:
    """Braille spinner with Claude Code-style colored output."""

    _FRAMES = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]

    def __init__(self, message: str = "Thinking"):
        self._message = message
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None

    def start(self) -> None:
        self._stop.clear()
        self._thread = threading.Thread(target=self._spin, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread:
            self._thread.join()
        clear_len = len(self._message) + 8
        sys.stdout.write("\r" + " " * clear_len + "\r")
        sys.stdout.flush()

    def _spin(self) -> None:
        idx = 0
        while not self._stop.is_set():
            frame = self._FRAMES[idx % len(self._FRAMES)]
            sys.stdout.write(
                f"\r  {C.YELLOW}{frame}{C.R} {C.DIM}{self._message}…{C.R}"
            )
            sys.stdout.flush()
            time.sleep(0.1)
            idx += 1


# ---------------------------------------------------------------------------
# Response display
# ---------------------------------------------------------------------------

def _label(icon: str, text: str) -> str:
    """Render a bold-cyan labeled header line: '◆ Text'."""
    return f"{C.LABEL}{icon} {text}{C.R}"


def _print_response(response: "ChemResponse", verbose: bool = False) -> None:
    """Print a ChemResponse in Claude Code-style formatted output."""
    import json

    print()

    # ── Plan / Hypothesis ──────────────────────────────────────────────
    if response.hypothesis or response.plan_rationale:
        conf_str = f"{C.CONF}[{response.confidence:.0%}]{C.R}" if response.confidence else ""
        print(_label("◆", "Hypothesis") + f"  {conf_str}")
        if response.hypothesis:
            print(f"  {C.HYPO}{response.hypothesis}{C.R}")
        if response.plan_rationale:
            print(f"  {C.META}{response.plan_rationale}{C.R}")
        if response.plan_revised:
            print(f"  {C.META}↺ Plan revised after Group 0 observation{C.R}")
        print()

    # ── Tools called ──────────────────────────────────────────────────
    if response.tools_called:
        # Arrow-joined tool chain
        arrow = f"  {C.META}→{C.R}  "
        tool_chain = arrow.join(f"{C.TOOL}{t}{C.R}" for t in response.tools_called)
        print(_label("⎿", "Tools"))
        print(f"  {tool_chain}")

        if verbose and response.tool_results:
            print()
            for name, result in response.tool_results.items():
                print(f"  {C.DIM}┌ {name}{C.R}")
                if isinstance(result, dict):
                    display = {k: v for k, v in result.items() if k != "success"}
                    try:
                        raw = json.dumps(display, indent=2, default=str)[:1000]
                    except Exception:
                        raw = str(display)[:600]
                    for line in raw.splitlines():
                        print(f"  {C.DIM}│{C.R}  {line}")
                else:
                    print(f"  {C.DIM}│{C.R}  {str(result)[:400]}")
    else:
        print(f"  {C.META}⎿ (no tools called — answered from LLM knowledge){C.R}")

    # ── Answer ────────────────────────────────────────────────────────
    print()
    print(_SEP_FAT)
    for line in response.answer.strip().splitlines():
        # Preserve blank lines; indent content lines
        print(f"  {line}" if line.strip() else "")
    print(_SEP_FAT)

    # ── Metadata footer ───────────────────────────────────────────────
    parts = [
        response.model,
        f"{response.elapsed_s:.1f}s",
        f"{response.llm_calls} LLM calls",
        f"{len(response.tools_called)} tools",
    ]
    print(f"  {C.META}{' · '.join(parts)}{C.R}")

    # ── Warnings ──────────────────────────────────────────────────────
    if response.warnings:
        print()
        for w in response.warnings:
            print(f"  {C.WARN}⚠  {w}{C.R}")


# ---------------------------------------------------------------------------
# Main CLI loop
# ---------------------------------------------------------------------------

def _cmd_setup(args: argparse.Namespace) -> None:  # noqa: ARG001
    """Handle the 'setup' subcommand — choose and save default model."""
    print()
    print(f"  {C.LABEL}◆ ChemCoworker Setup{C.R}  {C.DIM}Choose your default model{C.R}")
    print(f"  {_SEP}")

    # Show current saved model (if any)
    current = _load_config()
    current_label = f"{current['name']}  ({current['provider']})"
    print(f"  {C.META}Current:{C.R}  {C.DIM}{current_label}{C.R}")
    print()

    selected = select_model_interactive()
    _save_config(selected["name"], selected["provider"])

    provider_color = C.CYAN if selected["provider"] == "openai" else C.MAGENTA
    print()
    print(
        f"  {C.OK}✓{C.R}  Saved  "
        f"{C.BOLD}{selected['name']}{C.R}  "
        f"{provider_color}{selected['provider']}{C.R}"
    )
    print(f"  {C.META}Config:{C.R}  {C.DIM}{_CONFIG_PATH}{C.R}")
    print()


def _cmd_intake(args: argparse.Namespace) -> None:
    """Handle the 'intake' subcommand."""
    import time as _time
    from chem_coworker.extractor import NotesExtractor

    if args.model is not None:
        model = args.model
        provider = args.provider or ("aliyun" if model in _ALIYUN_MODELS else "openai")
    else:
        saved = _load_config()
        model = saved["name"]
        provider = args.provider or saved["provider"]

    print()
    print(f"  {C.LABEL}◆ ChemCoworker Intake{C.R}  {C.DIM}Document → Notes{C.R}")
    print(f"  {_SEP}")
    print(f"  {C.META}Model:{C.R}  {C.BOLD}{model}{C.R}  {C.META}{provider}{C.R}")
    print(f"  {C.META}Source:{C.R} {C.DIM}{args.source!r}{C.R}")
    if args.reaction_type:
        print(f"  {C.META}Hint:{C.R}   {C.DIM}{args.reaction_type}{C.R}")
    print()

    spinner = Spinner("Extracting notes")
    spinner.start()
    t0 = _time.time()

    try:
        extractor = NotesExtractor(provider=provider, model=model, verbose=args.verbose)
        result = extractor.intake(
            source=args.source,
            reaction_type=args.reaction_type or "",
            save_to_literature=not args.no_save,
        )
    except Exception as exc:
        spinner.stop()
        print(f"  {C.ERR}✗{C.R}  {exc}")
        sys.exit(1)

    spinner.stop()
    elapsed = _time.time() - t0

    if not result.get("success"):
        print(f"  {C.ERR}✗{C.R}  {result.get('error', 'Unknown error')}")
        for w in result.get("warnings", []):
            print(f"  {C.WARN}⚠{C.R}  {w}")
        sys.exit(1)

    # ── Results ──────────────────────────────────────────────────────
    print(f"  {C.OK}✓{C.R}  Extracted {result['char_count']:,} chars  {C.META}{elapsed:.1f}s{C.R}")
    print()

    # Reaction types detected
    types_str = ", ".join(result["reaction_types"])
    print(f"  {C.LABEL}◆ Reaction types{C.R}  {C.TOOL}{types_str}{C.R}")

    # Notes files written
    print(f"  {C.LABEL}◆ Notes written to{C.R}")
    for nf in result["notes_files"]:
        print(f"      {C.DIM}{nf}{C.R}")

    # Warnings
    for w in result.get("warnings", []):
        print(f"  {C.WARN}⚠{C.R}  {w}")

    # Preview extracted notes
    if not args.quiet:
        print()
        print(_SEP_FAT)
        preview = result["extracted_notes"]
        if len(preview) > 2000:
            preview = preview[:2000] + f"\n\n{C.DIM}[… {len(result['extracted_notes']) - 2000} more chars …]{C.R}"
        for line in preview.splitlines():
            print(f"  {line}" if line.strip() else "")
        print(_SEP_FAT)


def main() -> None:
    arg_parser = argparse.ArgumentParser(
        description="ChemCoworker — general-purpose chemistry AI agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python chem_coworker/cli.py
  python chem_coworker/cli.py --model gpt-5.2
  python chem_coworker/cli.py --model glm-5 --provider aliyun
  python chem_coworker/cli.py --verbose
  python chem_coworker/cli.py intake https://www.orgsyn.org/demo.aspx?prep=v102p0086
  python chem_coworker/cli.py intake paper.pdf --reaction-type buchwald_hartwig
        """,
    )
    arg_parser.add_argument("--model", default=None, help="LLM model name (skip selector if set)")
    arg_parser.add_argument("--provider", default=None, help="LLM provider: openai or aliyun")
    arg_parser.add_argument("--verbose", action="store_true", help="Show tool result details")

    subparsers = arg_parser.add_subparsers(dest="command")

    # ── setup subcommand ──────────────────────────────────────────────
    subparsers.add_parser(
        "setup",
        help="Choose and save your default model (written to ~/.chemcoworker/config.json)",
    )

    # ── intake subcommand ─────────────────────────────────────────────
    intake_parser = subparsers.add_parser(
        "intake",
        help="Extract chemistry notes from a URL, PDF, or text file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python chem_coworker/cli.py intake https://www.orgsyn.org/demo.aspx?prep=v102p0086
  python chem_coworker/cli.py intake paper.pdf --reaction-type suzuki_miyaura
  python chem_coworker/cli.py intake my_notes.txt
        """,
    )
    intake_parser.add_argument("source", help="URL, file path, or raw text to process")
    intake_parser.add_argument(
        "--reaction-type", default="",
        help="Hint for reaction type (e.g. suzuki_miyaura). Auto-detected if omitted.",
    )
    intake_parser.add_argument(
        "--no-save", action="store_true",
        help="Do not save fetched document to literature/ folder",
    )
    intake_parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress extracted notes preview",
    )

    args = arg_parser.parse_args()

    # ── Route to subcommands ──────────────────────────────────────────
    if args.command == "setup":
        _cmd_setup(args)
        return
    if args.command == "intake":
        _cmd_intake(args)
        return

    # ── Banner ────────────────────────────────────────────────────────
    print()
    print(f"  {C.LABEL}◆ ChemCoworker{C.R}  {C.DIM}Chemistry AI Agent{C.R}")
    print(f"  {C.META}General-purpose chemistry Q&A, analysis, and prediction{C.R}")
    print(f"  {_SEP}")

    # ── Model resolution (CLI flag > saved config > built-in default) ─
    if args.model is not None:
        model = args.model
        provider = args.provider or ("aliyun" if model in _ALIYUN_MODELS else "openai")
    else:
        saved = _load_config()
        model = saved["name"]
        provider = args.provider or saved["provider"]

    # ── API key check ─────────────────────────────────────────────────
    key_env = f"{provider.upper()}_API_KEY"
    if not os.getenv(key_env):
        print(f"\n  {C.ERR}✗{C.R}  {key_env} not set in environment.")
        print(f"     Set it with: export {key_env}=your_api_key")
        sys.exit(1)

    provider_color = C.CYAN if provider == "openai" else C.MAGENTA
    print(f"\n  {C.META}Using{C.R}  {C.BOLD}{model}{C.R}  {provider_color}{provider}{C.R}")
    print(f"  {C.META}Type{C.R}  {C.DIM}exit · clear · tools  (run 'setup' to change model){C.R}")

    # ── Example queries ───────────────────────────────────────────────
    print(f"\n  {C.META}Examples:{C.R}")
    examples = [
        "Recommend conditions: Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "Explain why Pd(0) is needed for Suzuki coupling",
        "What bases work for Buchwald-Hartwig C-N coupling?",
        "What are the properties of c1cccnc1? Is it drug-like?",
        "My Suzuki gave only 30% yield — what could cause this?",
    ]
    for ex in examples:
        print(f"  {C.DIM}•{C.R}  {C.DIM}{ex}{C.R}")
    print()

    # ── Initialize agent ──────────────────────────────────────────────
    try:
        from chem_coworker import ChemCoworker, REGISTRY
        coworker = ChemCoworker(provider=provider, model=model, verbose=args.verbose)
        print(f"  {C.OK}✓{C.R}  Agent ready  {C.META}{len(REGISTRY)} tools registered{C.R}")
    except Exception as exc:
        print(f"  {C.ERR}✗{C.R}  Failed to initialize: {exc}")
        sys.exit(1)

    # ── Conversation loop ─────────────────────────────────────────────
    history: List[Dict] = []

    while True:
        print()
        try:
            query = input(f"  {C.PROMPT}>{C.R} ").strip()
        except (KeyboardInterrupt, EOFError):
            print(f"\n\n  {C.META}Goodbye!{C.R}")
            break

        if not query:
            continue

        if query.lower() in ("exit", "quit", "q"):
            print(f"  {C.META}Goodbye!{C.R}")
            break

        if query.lower() == "clear":
            history = []
            print(f"  {C.META}History cleared.{C.R}")
            continue

        if query.lower() == "tools":
            print()
            print(REGISTRY.describe_tools())
            continue

        # ── Run the agent ─────────────────────────────────────────────
        spinner = Spinner("ChemCoworker thinking")
        spinner.start()
        try:
            response, history = coworker.chat(query, history)
        except Exception as exc:
            spinner.stop()
            print(f"\n  {C.ERR}✗{C.R}  {exc}")
            continue
        spinner.stop()

        _print_response(response, verbose=args.verbose)


if __name__ == "__main__":
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent.parent))
    main()
