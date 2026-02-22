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
    python chem_coworker/cli.py intake paper.pdf --extract-model deepseek-v3.2 --extract-provider aliyun

Features:
  - Saved model preference (~/.chemcoworker/config.json)
  - Multi-turn conversation with history
  - Shows plan, tools, and answer in Claude Code-style output
  - Spinner animation during LLM calls
  - Wall-clock timing display
  - Intake subcommand: extract chemistry notes from URLs, PDFs, and text files
"""
from __future__ import annotations

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated", category=UserWarning)

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

# When run directly (python chem_coworker/cli.py) the project root is not on
# sys.path. Insert it so that `import chem_coworker` resolves correctly.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from chem_coworker.event_bus import EventBus, ChemEvent


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

    if response.streamed:
        # A5 — answer was already printed token-by-token.
        # Flush any partial line, close the answer block, then print metadata.
        _stream_state["writer"].flush_remaining()
        print()          # ensure we're on a fresh line
        print(_SEP_FAT)  # closing separator

        if verbose and response.tool_results:
            # Show tool result details below the streamed answer
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
        # Non-streaming: print everything at once (original path)
        print()

        # ── Plan / Hypothesis ──────────────────────────────────────────────
        if response.hypothesis or response.plan_rationale:
            conf_str = f"{C.CONF}[{response.confidence:.0%}]{C.R}" if response.confidence else ""
            print(_label("◆", "Hypothesis") + f"  {conf_str}")
            if response.hypothesis:
                print(f"  {C.HYPO}{response.hypothesis}{C.R}")
            if response.plan_rationale:
                print(f"  {C.META}{response.plan_rationale}{C.R}")
            print()

        # ── Tools called ──────────────────────────────────────────────────
        if response.tools_called:
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
            print(f"  {line}" if line.strip() else "")
        print(_SEP_FAT)

    # ── Metadata footer (always printed) ─────────────────────────────
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
# A3 — Real-time tool streaming callbacks
# ---------------------------------------------------------------------------

_running_tools: set = set()       # tools currently executing
_phase_spinner: list = [None]     # mutable container for the active phase spinner


def _progress(event: str, name: str, elapsed: float) -> None:
    """Progress callback: prints one line per tool as it starts/completes."""
    if event == "start":
        _running_tools.add(name)
        sys.stdout.write(
            f"\r  {C.YELLOW}⠋{C.R} {C.DIM}{' · '.join(sorted(_running_tools))}…{C.R}   "
        )
        sys.stdout.flush()
    else:
        _running_tools.discard(name)
        icon = f"{C.OK}✓" if event == "done" else f"{C.ERR}✗"
        # Clear spinner line then print completed tool
        print(f"\r  {icon}{C.R}  {C.TOOL}{name}{C.R}  {C.META}{elapsed:.1f}s{C.R}     ")


_PHASE_LABELS = {
    "reason_start":  "Reasoning",
    "observe_start": "Observing",
    "synth_start":   "Synthesizing",
    "compact_start": "Compacting history",
}


# ---------------------------------------------------------------------------
# A5 — Answer streaming
# ---------------------------------------------------------------------------

class _LineWriter:
    """
    Buffers streaming tokens and writes them to stdout with a 2-space indent
    on each line — matching the non-streaming answer display style.
    """
    def __init__(self) -> None:
        self._buf = ""

    def write(self, token: str) -> None:
        self._buf += token
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            sys.stdout.write(f"  {line}\n")
        sys.stdout.flush()

    def flush_remaining(self) -> None:
        """Flush any partial line at end of stream (no trailing newline)."""
        if self._buf:
            sys.stdout.write(f"  {self._buf}")
            sys.stdout.flush()
            self._buf = ""

    def reset(self) -> None:
        self._buf = ""


# Mutable streaming state — reset before each query
_stream_state: dict = {
    "first_chunk": True,      # True until first token arrives
    "pre_synth_info": {},     # set by _pre_synth_cb before synthesis
    "writer": _LineWriter(),  # line-buffered stdout writer
}


def _print_pre_answer_chrome(ctx: dict) -> None:
    """Print hypothesis + tools block. Used both in streaming and non-streaming paths."""
    hypothesis  = ctx.get("hypothesis", "")
    confidence  = ctx.get("confidence", 0.0)
    rationale   = ctx.get("rationale", "")
    tools_called = ctx.get("tools_called", [])

    print()
    if hypothesis or rationale:
        conf_str = f"{C.CONF}[{confidence:.0%}]{C.R}" if confidence else ""
        print(_label("◆", "Hypothesis") + f"  {conf_str}")
        if hypothesis:
            print(f"  {C.HYPO}{hypothesis}{C.R}")
        if rationale:
            print(f"  {C.META}{rationale}{C.R}")
        print()

    if tools_called:
        arrow = f"  {C.META}→{C.R}  "
        tool_chain = arrow.join(f"{C.TOOL}{t}{C.R}" for t in tools_called)
        print(_label("⎿", "Tools"))
        print(f"  {tool_chain}")
    else:
        print(f"  {C.META}⎿ (no tools called — answered from LLM knowledge){C.R}")
    print()


def _pre_synth_cb(
    hypothesis: str,
    confidence: float,
    rationale: str,
    tools_called: list,
) -> None:
    """Called by agent right before synthesis starts — stores plan info for streaming."""
    _stream_state["pre_synth_info"] = {
        "hypothesis": hypothesis,
        "confidence": confidence,
        "rationale": rationale,
        "tools_called": tools_called,
    }


def _stream_token(token: str) -> None:
    """Called by agent for each synthesis token. Prints pre-answer chrome on first call."""
    if _stream_state["first_chunk"]:
        _stream_state["first_chunk"] = False
        # Stop the Synthesizing... spinner before printing anything
        if _phase_spinner[0] is not None:
            _phase_spinner[0].stop()
            _phase_spinner[0] = None
        # Print hypothesis + tools block
        _print_pre_answer_chrome(_stream_state["pre_synth_info"])
        # Opening answer separator
        print(_SEP_FAT)
        _stream_state["writer"].reset()

    _stream_state["writer"].write(token)


def _phase(phase: str) -> None:
    """Phase callback: starts/stops spinners around LLM calls."""
    if phase.endswith("_start") and phase in _PHASE_LABELS:
        _phase_spinner[0] = Spinner(_PHASE_LABELS[phase])
        _phase_spinner[0].start()
    elif phase.endswith("_done") and _phase_spinner[0] is not None:
        _phase_spinner[0].stop()
        _phase_spinner[0] = None


# ---------------------------------------------------------------------------
# A2 — Plan approval helpers
# ---------------------------------------------------------------------------

def _show_plan_and_confirm(plan):
    """
    Display the structured plan and pause for user confirmation (--plan mode).
    Returns the (possibly modified) plan or raises PlanRejected.
    """
    from chem_coworker.plan import PlanRejected

    print()
    conf_str = f"{C.CONF}[{plan.confidence:.0%}]{C.R}" if plan.confidence else ""
    print(_label("◆", "Proposed Plan") + f"  {conf_str}")
    if plan.hypothesis:
        print(f"  {C.HYPO}{plan.hypothesis}{C.R}")
    if plan.rationale:
        print(f"  {C.META}{plan.rationale}{C.R}")
    print()
    for i, group in enumerate(plan.groups):
        names = "  ·  ".join(f"{C.TOOL}{c.name}{C.R}" for c in group)
        print(f"  {C.DIM}Group {i}{C.R}  {names}")
    print()

    try:
        raw = input(f"  {C.PROMPT}>{C.R} Proceed? [Y/n/skip <tool>]: ").strip().lower()
    except (KeyboardInterrupt, EOFError):
        raise PlanRejected("Cancelled by user.")

    if raw in ("n", "no", "q", "quit"):
        raise PlanRejected("Cancelled by user.")
    if raw.startswith("skip "):
        tool_name = raw[5:].strip()
        return _drop_tool_from_plan(plan, tool_name)
    print()  # blank line before streaming starts
    return plan


def _drop_tool_from_plan(plan, name: str):
    """Return a copy of plan with the named tool removed from all groups."""
    import dataclasses
    new_groups = [
        [c for c in group if c.name != name]
        for group in plan.groups
    ]
    new_groups = [g for g in new_groups if g]
    return dataclasses.replace(plan, groups=new_groups)


# ---------------------------------------------------------------------------
# In-session coworker factory
# ---------------------------------------------------------------------------

def _build_event_bus() -> EventBus:
    """Wire CLI display functions to a fresh EventBus for one session."""
    bus = EventBus()

    # Tool progress — TOOL_START shows spinner; TOOL_DONE / TOOL_ERROR prints result
    bus.subscribe(ChemEvent.TOOL_START, lambda tool_name, **_: _progress("start", tool_name, 0.0))
    bus.subscribe(ChemEvent.TOOL_DONE, lambda tool_name, elapsed_s, **_: _progress("done", tool_name, elapsed_s))
    bus.subscribe(ChemEvent.TOOL_ERROR, lambda tool_name, elapsed_s, **_: _progress("error", tool_name, elapsed_s))

    # Phase spinners — reuse the existing _PHASE_LABELS mapping
    bus.subscribe(ChemEvent.PHASE_START, lambda phase, **_: _phase(f"{phase}_start"))
    bus.subscribe(ChemEvent.PHASE_END, lambda phase, **_: _phase(f"{phase}_done"))

    # Compact spinners
    bus.subscribe(ChemEvent.COMPACT_START, lambda **_: _phase("compact_start"))
    bus.subscribe(ChemEvent.COMPACT_END, lambda **_: _phase("compact_done"))

    # Answer streaming
    bus.subscribe(ChemEvent.PRE_SYNTH, lambda hypothesis, confidence, rationale, tools_called, **_:
        _pre_synth_cb(hypothesis, confidence, rationale, tools_called))
    bus.subscribe(ChemEvent.STREAM_TOKEN, lambda token, **_: _stream_token(token))

    return bus


def _init_coworker(model: str, provider: str, verbose: bool, plan_mode: bool):
    """Create a ChemCoworker with the current session settings."""
    from chem_coworker import ChemCoworker
    return ChemCoworker(
        provider=provider,
        model=model,
        verbose=verbose,
        event_bus=_build_event_bus(),
        plan_callback=_show_plan_and_confirm if plan_mode else None,
    )


def _print_settings(model: str, provider: str, verbose: bool, plan_mode: bool) -> None:
    """Print current session settings."""
    provider_color = C.CYAN if provider == "openai" else C.MAGENTA
    print()
    print(_label("◆", "Current Settings"))
    print(f"  {C.META}Model   {C.R}  {C.BOLD}{model}{C.R}  {provider_color}{provider}{C.R}")
    print(f"  {C.META}Verbose {C.R}  {C.OK if verbose else C.DIM}{'on' if verbose else 'off'}{C.R}")
    print(f"  {C.META}Plan    {C.R}  {C.OK if plan_mode else C.DIM}{'on' if plan_mode else 'off'}{C.R}")
    print(f"  {C.META}Config  {C.R}  {C.DIM}{_CONFIG_PATH}{C.R}")
    print()
    print(f"  {C.DIM}Commands: /model · /plan · /verbose · /compact · /settings{C.R}")


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

    # Main model (fallback when --extract-model not set)
    if args.model is not None:
        model = args.model
    else:
        saved = _load_config()
        model = saved["name"]

    # Extraction model — can be different (cheaper/faster) than the chat model
    extract_model = args.extract_model or model
    extract_provider = (
        args.extract_provider
        or ("aliyun" if extract_model in _ALIYUN_MODELS else "openai")
    )

    print()
    print(f"  {C.LABEL}◆ ChemCoworker Intake{C.R}  {C.DIM}Document → Notes{C.R}")
    print(f"  {_SEP}")
    if extract_model != model:
        print(f"  {C.META}Extract model:{C.R}  {C.BOLD}{extract_model}{C.R}  {C.META}{extract_provider}{C.R}")
    else:
        print(f"  {C.META}Model:{C.R}  {C.BOLD}{extract_model}{C.R}  {C.META}{extract_provider}{C.R}")
    print(f"  {C.META}Source:{C.R} {C.DIM}{args.source!r}{C.R}")
    if args.reaction_type:
        print(f"  {C.META}Hint:{C.R}   {C.DIM}{args.reaction_type}{C.R}")
    note_type = getattr(args, "note_type", "reactions") or "reactions"
    print(f"  {C.META}Type:{C.R}   {C.DIM}{note_type}{C.R}")
    print()

    spinner = Spinner("Extracting notes")
    spinner.start()
    t0 = _time.time()

    try:
        extractor = NotesExtractor(provider=extract_provider, model=extract_model, verbose=args.verbose)
        result = extractor.intake(
            source=args.source,
            reaction_type=args.reaction_type or "",
            note_type=note_type,
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
    arg_parser.add_argument(
        "--plan", action="store_true",
        help="A2: Show proposed tool plan and confirm before executing (plan approval mode)",
    )

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
        "--note-type", default="reactions",
        choices=["reactions", "mechanisms", "substrates", "protocols", "routes"],
        help="Which notes subfolder to write to (default: reactions).",
    )
    intake_parser.add_argument(
        "--extract-model", default=None,
        help=(
            "Model to use for extraction (overrides --model for this step). "
            "Useful for a cheaper/faster model: e.g. --extract-model deepseek-v3.2 "
            "--extract-provider aliyun"
        ),
    )
    intake_parser.add_argument(
        "--extract-provider", default=None,
        help="Provider for --extract-model (openai or aliyun). Auto-detected if omitted.",
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

    # Mutable session state (can be changed mid-session via commands)
    verbose   = args.verbose
    plan_mode = args.plan

    provider_color = C.CYAN if provider == "openai" else C.MAGENTA
    print(f"\n  {C.META}Using{C.R}  {C.BOLD}{model}{C.R}  {provider_color}{provider}{C.R}")
    print(f"  {C.META}Type{C.R}  {C.DIM}/model · /plan · /verbose · /settings · exit{C.R}")

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
    from chem_coworker import REGISTRY
    try:
        coworker = _init_coworker(model, provider, verbose, plan_mode)
        print(f"  {C.OK}✓{C.R}  Agent ready  {C.META}{len(REGISTRY)} tools registered{C.R}")
        if plan_mode:
            print(f"  {C.META}Plan approval ON{C.R}")
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

        ql = query.lower()

        # ── Built-in commands ─────────────────────────────────────────
        if ql in ("exit", "quit", "q"):
            print(f"  {C.META}Goodbye!{C.R}")
            break

        if ql == "clear":
            history = []
            print(f"  {C.META}History cleared.{C.R}")
            continue

        if ql == "tools":
            print()
            print(REGISTRY.describe_tools())
            continue

        # ── Settings commands (/model, /plan, /verbose, /settings, /compact)
        if ql in ("/model", "change model", "model"):
            selected = select_model_interactive()
            new_model    = selected["name"]
            new_provider = selected["provider"]
            _save_config(new_model, new_provider)
            try:
                coworker = _init_coworker(new_model, new_provider, verbose, plan_mode)
                model, provider = new_model, new_provider
                pc = C.CYAN if provider == "openai" else C.MAGENTA
                print(f"  {C.OK}✓{C.R}  Switched to {C.BOLD}{model}{C.R}  {pc}{provider}{C.R}")
                print(f"  {C.META}Config saved. History preserved.{C.R}")
            except Exception as exc:
                print(f"  {C.ERR}✗{C.R}  Failed to switch model: {exc}")
            continue

        if ql in ("/plan", "toggle plan"):
            plan_mode = not plan_mode
            coworker.plan_callback = _show_plan_and_confirm if plan_mode else None
            state = f"{C.OK}ON{C.R}" if plan_mode else f"{C.DIM}off{C.R}"
            print(f"  {C.META}Plan approval{C.R}  {state}")
            continue

        if ql in ("/verbose", "toggle verbose"):
            verbose = not verbose
            coworker.verbose = verbose
            coworker.executor.verbose = verbose
            state = f"{C.OK}ON{C.R}" if verbose else f"{C.DIM}off{C.R}"
            print(f"  {C.META}Verbose{C.R}  {state}")
            continue

        if ql in ("/settings", "settings"):
            _print_settings(model, provider, verbose, plan_mode)
            continue

        if ql in ("/compact", "compact"):
            if len(history) < 2:
                print(f"  {C.META}Nothing to compact (history is empty).{C.R}")
            else:
                history = coworker._compact_history(history)
                print(f"  {C.META}↺ History compacted to {len(history)} messages.{C.R}")
            continue

        # ── Run the agent ─────────────────────────────────────────────
        # A3/A5: clear streaming state before each query
        _running_tools.clear()
        _stream_state["first_chunk"] = True
        _stream_state["pre_synth_info"] = {}
        _stream_state["writer"].reset()
        try:
            response, history = coworker.chat(query, history)
        except Exception as exc:
            # Stop any lingering phase spinner on error
            if _phase_spinner[0] is not None:
                _phase_spinner[0].stop()
                _phase_spinner[0] = None
            print(f"\n  {C.ERR}✗{C.R}  {exc}")
            continue

        # A4: show compaction notice if history was auto-summarized this turn
        if response.compacted:
            print(f"  {C.META}↺ Conversation history compacted{C.R}")

        _print_response(response, verbose=verbose)


if __name__ == "__main__":
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent.parent))
    main()
