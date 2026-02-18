"""
ChemCoworker Interactive CLI

Usage:
    python chem_coworker/cli.py
    python chem_coworker/cli.py --model gpt-5.2 --provider openai
    python chem_coworker/cli.py --model glm-5 --provider aliyun

Features:
  - Interactive model selector at startup
  - Multi-turn conversation with history
  - Shows [PLAN], [TOOLS CALLED], and [ANSWER] blocks
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
# Model selector (same pattern as existing CLIs)
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
    print("\n┌─ Select LLM Model " + "─" * 40)
    for i, m in enumerate(SELECTABLE_MODELS, 1):
        tag = " ← default" if i == 1 else ""
        print(f"│  {i}. {m['name']:<20} [{m['provider']}]{tag}")
    print("└" + "─" * 59)
    raw = input("  Choice (1–{}, Enter=default): ".format(len(SELECTABLE_MODELS))).strip()
    if not raw:
        return SELECTABLE_MODELS[0]
    try:
        idx = int(raw) - 1
        if 0 <= idx < len(SELECTABLE_MODELS):
            return SELECTABLE_MODELS[idx]
    except ValueError:
        pass
    print("  Invalid selection — using default (o4-mini)")
    return SELECTABLE_MODELS[0]


# ---------------------------------------------------------------------------
# Spinner
# ---------------------------------------------------------------------------

class Spinner:
    """Simple terminal spinner for async feedback."""

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
        sys.stdout.write("\r" + " " * (len(self._message) + 4) + "\r")
        sys.stdout.flush()

    def _spin(self) -> None:
        idx = 0
        while not self._stop.is_set():
            frame = self._FRAMES[idx % len(self._FRAMES)]
            sys.stdout.write(f"\r{frame} {self._message}...")
            sys.stdout.flush()
            time.sleep(0.1)
            idx += 1


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

_SEP = "─" * 60


def _print_section(title: str, content: str, color_code: str = "") -> None:
    reset = "\033[0m" if color_code else ""
    print(f"\n{color_code}┌─ {title} {'─' * max(0, 55 - len(title))}{reset}")
    for line in content.strip().splitlines():
        print(f"{color_code}│{reset}  {line}")
    print(f"{color_code}└{reset}" + _SEP)


def _print_response(response: "ChemResponse", verbose: bool = False) -> None:
    """Print a ChemResponse in a structured, readable format."""
    # Timing line
    timing = f"  [{response.model} | {response.elapsed_s:.1f}s | {response.llm_calls} LLM calls | {len(response.tools_called)} tools]"

    # [PLAN] section — show hypothesis + rationale (not full plan text)
    if response.hypothesis or response.plan_rationale:
        plan_display = ""
        if response.hypothesis:
            plan_display += f"Hypothesis: {response.hypothesis}\n"
        if response.plan_rationale:
            plan_display += f"Rationale:  {response.plan_rationale}"
        _print_section("PLAN", plan_display, "\033[36m")  # cyan

    # [TOOLS CALLED] section
    if response.tools_called:
        tools_display = " → ".join(response.tools_called)
        _print_section("TOOLS CALLED", tools_display, "\033[33m")  # yellow
        if verbose and response.tool_results:
            import json
            for name, result in response.tool_results.items():
                print(f"    [{name}]:")
                if isinstance(result, dict):
                    display = {k: v for k, v in result.items() if k not in ("success",)}
                    try:
                        print("    " + json.dumps(display, indent=2, default=str)[:800].replace("\n", "\n    "))
                    except Exception:
                        print(f"    {str(result)[:400]}")
    else:
        print("\n\033[33m│\033[0m  (No tools called — answered from LLM knowledge)")

    # [ANSWER] section
    _print_section("ANSWER", response.answer, "\033[32m")  # green
    print(timing)

    # Warnings
    if response.warnings:
        print("\n  \033[31m⚠ Warnings:\033[0m")
        for w in response.warnings:
            print(f"    • {w}")


# ---------------------------------------------------------------------------
# Main CLI loop
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
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
    parser.add_argument("--model", default=None, help="LLM model name (skip selector if set)")
    parser.add_argument("--provider", default=None, help="LLM provider: openai or aliyun")
    parser.add_argument("--verbose", action="store_true", help="Show tool result details")
    args = parser.parse_args()

    # ── Banner ──────────────────────────────────────────────────────────
    print("\n" + "═" * 60)
    print("  ChemCoworker — Chemistry AI Agent")
    print("  General-purpose chemistry Q&A, analysis, and prediction")
    print("=" * 60)

    # ── Model selection ─────────────────────────────────────────────────
    if args.model is None:
        selected = select_model_interactive()
        model = selected["name"]
        provider = selected["provider"]
    else:
        model = args.model
        provider = args.provider or ("aliyun" if model in _ALIYUN_MODELS else "openai")

    # ── API key check ────────────────────────────────────────────────────
    key_env = f"{provider.upper()}_API_KEY"
    if not os.getenv(key_env):
        print(f"\n  ✗ {key_env} not set in environment.")
        print(f"  Set it with: export {key_env}=your_api_key")
        sys.exit(1)

    print(f"\n  Using: {model} [{provider}]")
    print("  Type 'exit' or 'quit' to leave | 'clear' to reset history")
    print("  Type 'tools' to list available tools")
    print(_SEP)

    # ── Example queries ──────────────────────────────────────────────────
    print("\n  Example queries:")
    print("    • Recommend conditions: Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
    print("    • Explain why Pd(0) is needed for Suzuki coupling")
    print("    • What bases work for Buchwald-Hartwig C-N coupling?")
    print("    • What are the properties of c1cccnc1? Is it drug-like?")
    print("    • My Suzuki gave only 30% yield — what could cause this?")
    print()

    # ── Initialize agent ─────────────────────────────────────────────────
    try:
        from chem_coworker import ChemCoworker, REGISTRY
        coworker = ChemCoworker(provider=provider, model=model, verbose=args.verbose)
        print(f"  ✓ Agent ready | {len(REGISTRY)} tools registered")
    except Exception as exc:
        print(f"  ✗ Failed to initialize agent: {exc}")
        sys.exit(1)

    # ── Conversation loop ─────────────────────────────────────────────────
    history = []

    while True:
        print()
        try:
            query = input("You: ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\n\n  Goodbye!")
            break

        if not query:
            continue

        if query.lower() in ("exit", "quit", "q"):
            print("  Goodbye!")
            break

        if query.lower() == "clear":
            history = []
            print("  History cleared.")
            continue

        if query.lower() == "tools":
            print("\n" + REGISTRY.describe_tools())
            continue

        # ── Run the agent ────────────────────────────────────────────────
        spinner = Spinner("ChemCoworker thinking")
        spinner.start()
        try:
            response, history = coworker.chat(query, history)
        except Exception as exc:
            spinner.stop()
            print(f"\n  ✗ Error: {exc}")
            continue
        spinner.stop()

        _print_response(response, verbose=args.verbose)


if __name__ == "__main__":
    # Allow running as: python chem_coworker/cli.py
    # Add parent dir to path if running directly
    sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent.parent))
    main()
