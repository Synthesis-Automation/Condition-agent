"""
ChemCoworker Interactive CLI

Usage:
    python chem_coworker/cli.py
    python chem_coworker/cli.py --model gpt-5.2 --provider openai
    python chem_coworker/cli.py --model glm-5 --provider aliyun

Features:
  - Interactive model selector at startup
  - Multi-turn conversation with history
  - Shows plan, tools, and answer in Claude Code-style output
  - Spinner animation during LLM calls
  - Wall-clock timing display
"""
from __future__ import annotations

import argparse
import os
import sys
import threading
import time
from typing import Dict, List, Optional


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
        """,
    )
    arg_parser.add_argument("--model", default=None, help="LLM model name (skip selector if set)")
    arg_parser.add_argument("--provider", default=None, help="LLM provider: openai or aliyun")
    arg_parser.add_argument("--verbose", action="store_true", help="Show tool result details")
    args = arg_parser.parse_args()

    # ── Banner ────────────────────────────────────────────────────────
    print()
    print(f"  {C.LABEL}◆ ChemCoworker{C.R}  {C.DIM}Chemistry AI Agent{C.R}")
    print(f"  {C.META}General-purpose chemistry Q&A, analysis, and prediction{C.R}")
    print(f"  {_SEP}")

    # ── Model selection ───────────────────────────────────────────────
    if args.model is None:
        selected = select_model_interactive()
        model = selected["name"]
        provider = selected["provider"]
    else:
        model = args.model
        provider = args.provider or ("aliyun" if model in _ALIYUN_MODELS else "openai")

    # ── API key check ─────────────────────────────────────────────────
    key_env = f"{provider.upper()}_API_KEY"
    if not os.getenv(key_env):
        print(f"\n  {C.ERR}✗{C.R}  {key_env} not set in environment.")
        print(f"     Set it with: export {key_env}=your_api_key")
        sys.exit(1)

    provider_color = C.CYAN if provider == "openai" else C.MAGENTA
    print(f"\n  {C.META}Using{C.R}  {C.BOLD}{model}{C.R}  {provider_color}{provider}{C.R}")
    print(f"  {C.META}Type{C.R}  {C.DIM}exit · clear · tools{C.R}")

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
