"""
Main terminal CLI app (interactive REPL + subcommands).

How to use:
    # Recommended (package entrypoint)
    python -m chem_coworker._cli.app

    # Direct script execution (supported)
    python chem_coworker/_cli/app.py

    # Non-interactive single query
    python -m chem_coworker._cli.app ask "Explain why Pd(0) is needed for Suzuki coupling"

    # JSON output for automation
    python -m chem_coworker._cli.app ask "..." --output-format json

    # Batch mode (text lines or .jsonl)
    python -m chem_coworker._cli.app batch queries.txt --output-format jsonl
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from pathlib import Path
from typing import List, Optional

warnings.filterwarnings("ignore", message="pkg_resources is deprecated", category=UserWarning)

try:
    from .commands import (
        COMMAND_EXIT,
        COMMAND_HANDLED,
        build_default_command_registry,
        ReplSession,
    )
    from .config import CONFIG_PATH, load_config, save_config
    from .models import infer_provider, select_model_interactive
    from .renderers import render_ask_response, render_batch_failure, render_batch_success, render_batch_summary
    from .style import C, SEP
    from .ui import TerminalUI
except ImportError:
    # Support direct script execution: `python chem_coworker/_cli/app.py`
    _HERE = Path(__file__).resolve()
    _ROOT = _HERE.parents[2]
    if str(_ROOT) not in sys.path:
        sys.path.insert(0, str(_ROOT))
    from chem_coworker._cli.commands import (
        COMMAND_EXIT,
        COMMAND_HANDLED,
        build_default_command_registry,
        ReplSession,
    )
    from chem_coworker._cli.config import CONFIG_PATH, load_config, save_config
    from chem_coworker._cli.models import infer_provider, select_model_interactive
    from chem_coworker._cli.renderers import (
        render_ask_response,
        render_batch_failure,
        render_batch_success,
        render_batch_summary,
    )
    from chem_coworker._cli.style import C, SEP
    from chem_coworker._cli.ui import TerminalUI


def _init_coworker(
    model: str,
    provider: str,
    verbose: bool,
    plan_mode: bool,
    ui: Optional[TerminalUI] = None,
):
    from chem_coworker import ChemCoworker

    return ChemCoworker(
        provider=provider,
        model=model,
        verbose=verbose,
        event_bus=ui.build_event_bus() if ui is not None else None,
        plan_callback=ui.show_plan_and_confirm if (plan_mode and ui is not None) else None,
    )


def _resolve_model_provider(args: argparse.Namespace) -> tuple[str, str]:
    if args.model is not None:
        model = args.model
        provider = args.provider or infer_provider(model)
        return model, provider

    saved = load_config()
    model = saved["name"]
    provider = args.provider or saved["provider"]
    return model, provider


def _provider_color(provider: str) -> str:
    return C.CYAN if provider == "openai" else C.MAGENTA


def _print_api_key_hint(key_env: str) -> None:
    print(f"\n  {C.ERR}✗{C.R}  {key_env} not set in environment.")
    if os.name == "nt":
        print(f"     PowerShell: $env:{key_env}='your_api_key'")
    else:
        print(f"     Set it with: export {key_env}=your_api_key")


def _require_api_key(provider: str) -> None:
    key_env = f"{provider.upper()}_API_KEY"
    if not os.getenv(key_env):
        _print_api_key_hint(key_env)
        sys.exit(1)


def _load_batch_queries(path: str, input_format: str) -> List[str]:
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"Batch input file not found: {path}")

    text = file_path.read_text(encoding="utf-8")
    fmt = input_format
    if fmt == "auto":
        fmt = "jsonl" if file_path.suffix.lower() == ".jsonl" else "lines"

    queries: List[str] = []
    if fmt == "lines":
        for line in text.splitlines():
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            queries.append(raw)
    elif fmt == "jsonl":
        for idx, line in enumerate(text.splitlines(), 1):
            raw = line.strip()
            if not raw:
                continue
            try:
                item = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSONL at line {idx}: {exc}") from exc
            if isinstance(item, str):
                query = item.strip()
            elif isinstance(item, dict) and "query" in item:
                query = str(item["query"]).strip()
            else:
                raise ValueError(f"Invalid JSONL record at line {idx}: expected string or object with 'query'")
            if query:
                queries.append(query)
    else:
        raise ValueError(f"Unsupported batch format: {fmt}")

    return queries


def _cmd_ask(args: argparse.Namespace) -> None:
    from chem_coworker.response import ChemResponse

    query = " ".join(args.query).strip()
    if not query:
        print(f"  {C.ERR}✗{C.R}  Query is empty.")
        sys.exit(2)

    model, provider = _resolve_model_provider(args)
    _require_api_key(provider)

    output_format = args.output_format
    ui = None if output_format == "json" else TerminalUI()
    try:
        coworker = _init_coworker(model, provider, args.verbose, args.plan, ui)
        response: ChemResponse = coworker.run(query)
    except Exception as exc:
        print(f"  {C.ERR}✗{C.R}  {exc}")
        sys.exit(1)

    render_ask_response(response, output_format=output_format, verbose=args.verbose, ui=ui)


def _cmd_batch(args: argparse.Namespace) -> None:
    model, provider = _resolve_model_provider(args)
    _require_api_key(provider)
    output_format = args.output_format

    try:
        queries = _load_batch_queries(args.input, args.input_format)
    except Exception as exc:
        print(f"  {C.ERR}✗{C.R}  {exc}")
        sys.exit(1)

    if not queries:
        print(f"  {C.WARN}⚠{C.R}  No queries found in {args.input}")
        return
    if args.output and output_format != "jsonl":
        print(f"  {C.ERR}✗{C.R}  `--output` is supported only with `--output-format jsonl`.")
        sys.exit(2)

    ui = None if output_format == "jsonl" else TerminalUI()
    try:
        coworker = _init_coworker(model, provider, args.verbose, args.plan, ui)
    except Exception as exc:
        print(f"  {C.ERR}✗{C.R}  Failed to initialize: {exc}")
        sys.exit(1)

    out_fh = open(args.output, "w", encoding="utf-8") if args.output else sys.stdout
    close_output = bool(args.output)
    t0 = time.time()
    succeeded = 0
    failed = 0

    try:
        for idx, query in enumerate(queries, 1):
            try:
                response = coworker.run(query)
            except Exception as exc:
                failed += 1
                render_batch_failure(
                    idx=idx,
                    query=query,
                    error=exc,
                    output_format=output_format,
                    out_fh=out_fh,
                )
                if args.fail_fast or not args.continue_on_error:
                    break
                continue

            succeeded += 1
            render_batch_success(
                idx=idx,
                total=len(queries),
                query=query,
                response=response,
                output_format=output_format,
                verbose=args.verbose,
                ui=ui,
                out_fh=out_fh,
            )
    finally:
        if close_output:
            out_fh.close()

    elapsed = time.time() - t0
    render_batch_summary(
        total=len(queries),
        succeeded=succeeded,
        failed=failed,
        elapsed_s=elapsed,
        output_format=output_format,
    )
    if failed:
        sys.exit(1)


def _cmd_setup(args: argparse.Namespace) -> None:  # noqa: ARG001
    print()
    print(f"  {C.LABEL}◆ ChemCoworker Setup{C.R}  {C.DIM}Choose your default model{C.R}")
    print(f"  {SEP}")

    current = load_config()
    current_label = f"{current['name']}  ({current['provider']})"
    print(f"  {C.META}Current:{C.R}  {C.DIM}{current_label}{C.R}")
    print()

    selected = select_model_interactive()
    save_config(selected["name"], selected["provider"])

    provider_color = _provider_color(selected["provider"])
    print()
    print(
        f"  {C.OK}✓{C.R}  Saved  "
        f"{C.BOLD}{selected['name']}{C.R}  "
        f"{provider_color}{selected['provider']}{C.R}"
    )
    print(f"  {C.META}Config:{C.R}  {C.DIM}{CONFIG_PATH}{C.R}")
    print()


def _build_parser() -> argparse.ArgumentParser:
    arg_parser = argparse.ArgumentParser(
        description="ChemCoworker — general-purpose chemistry AI agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m chem_coworker._cli.app
  python -m chem_coworker._cli.app --model gpt-5.4
  python -m chem_coworker._cli.app --model glm-5 --provider aliyun
  python -m chem_coworker._cli.app --verbose
        """,
    )
    arg_parser.add_argument("--model", default=None, help="LLM model name (skip selector if set)")
    arg_parser.add_argument("--provider", default=None, help="LLM provider: openai or aliyun")
    arg_parser.add_argument("--verbose", action="store_true", help="Show tool result details")
    arg_parser.add_argument(
        "--plan",
        action="store_true",
        help="Show proposed tool plan and confirm before executing (plan approval mode)",
    )

    subparsers = arg_parser.add_subparsers(dest="command")
    subparsers.add_parser(
        "setup",
        help="Choose and save your default model (written to ~/.chemcoworker/config.json)",
    )

    ask_parser = subparsers.add_parser(
        "ask",
        help="Run a single non-interactive query",
    )
    ask_parser.add_argument("query", nargs="+", help="The question/query to run")
    ask_parser.add_argument(
        "--json",
        action="store_true",
        help="Deprecated alias for --output-format json",
    )
    ask_parser.add_argument(
        "--output-format",
        default="plain",
        choices=["plain", "json"],
        help="Output format (default: plain)",
    )

    batch_parser = subparsers.add_parser(
        "batch",
        help="Run multiple non-interactive queries from a file",
    )
    batch_parser.add_argument("input", help="Input file path (text lines or .jsonl)")
    batch_parser.add_argument(
        "--input-format",
        default="auto",
        choices=["auto", "lines", "jsonl"],
        help="Batch file format (default: auto by file extension)",
    )
    batch_parser.add_argument(
        "--json",
        action="store_true",
        help="Deprecated alias for --output-format jsonl",
    )
    batch_parser.add_argument(
        "--output-format",
        default="plain",
        choices=["plain", "jsonl"],
        help="Output format (default: plain)",
    )
    batch_parser.add_argument(
        "--output",
        default=None,
        help="Write batch output to file (stdout by default)",
    )
    batch_parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop processing on first failed query",
    )
    batch_parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue processing after errors (explicit CI-friendly mode; default behavior)",
    )
    return arg_parser


def _create_repl_session(
    coworker,
    model: str,
    provider: str,
    verbose: bool,
    plan_mode: bool,
    ui: TerminalUI,
    tool_registry,
):
    return ReplSession(
        coworker=coworker,
        model=model,
        provider=provider,
        verbose=verbose,
        plan_mode=plan_mode,
        condition_mode="auto",
        history=[],
        tool_registry=tool_registry,
        ui=ui,
        create_coworker=lambda m, p, v, pm: _init_coworker(m, p, v, pm, ui),
        save_default_config=save_config,
    )


def _is_condition_like_query(query: str) -> bool:
    text = str(query or "").lower()
    keywords = (
        "condition",
        "conditions",
        "catalyst",
        "ligand",
        "base",
        "solvent",
        "temperature",
        "reagent screen",
        "hte",
    )
    return any(k in text for k in keywords)


def _apply_condition_mode_hint(query: str, condition_mode: str) -> str:
    mode = str(condition_mode or "auto").strip().lower()
    if mode == "balanced":
        mode = "full"
    if mode not in {"auto", "full"}:
        mode = "auto"
    if mode == "auto" or not _is_condition_like_query(query):
        return query

    hint = (
        "[REPL condition mode preference: For condition recommendations, use "
        "condition_strategy='full' (cross-source analysis across literature, motif, "
        "similarity, and rules) unless the user explicitly asks for a single source.]"
    )
    return f"{query}\n\n{hint}"


def main(argv: Optional[List[str]] = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Backward-compatible aliases
    if getattr(args, "command", None) == "ask" and getattr(args, "json", False):
        args.output_format = "json"
    if getattr(args, "command", None) == "batch" and getattr(args, "json", False):
        args.output_format = "jsonl"
    if getattr(args, "command", None) == "batch":
        args.continue_on_error = bool(getattr(args, "continue_on_error", False) or not getattr(args, "fail_fast", False))

    if args.command == "setup":
        _cmd_setup(args)
        return
    if args.command == "ask":
        _cmd_ask(args)
        return
    if args.command == "batch":
        _cmd_batch(args)
        return
    print()
    print(f"  {C.LABEL}◆ ChemCoworker{C.R}  {C.DIM}Chemistry AI Agent{C.R}")
    print(f"  {C.META}General-purpose chemistry Q&A, analysis, and prediction{C.R}")
    print(f"  {SEP}")

    model, provider = _resolve_model_provider(args)
    _require_api_key(provider)

    verbose = args.verbose
    plan_mode = args.plan
    provider_color = _provider_color(provider)
    print(f"\n  {C.META}Using{C.R}  {C.BOLD}{model}{C.R}  {provider_color}{provider}{C.R}")
    print(f"  {C.META}Type{C.R}  {C.DIM}/help · /model · /plan · /verbose · /condmode · /settings · exit{C.R}")

    print(f"\n  {C.META}Examples:{C.R}")
    examples = [
        "Recommend conditions: Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "Recommend motif-screen conditions for my Suzuki coupling",
        "Find similar precedent conditions for Brc1ccc(F)cc1.OB(O)c1ccncc1>>...",
        "Explain why Pd(0) is needed for Suzuki coupling",
        "What bases work for Buchwald-Hartwig C-N coupling?",
        "What are the properties of c1cccnc1? Is it drug-like?",
        "My Suzuki gave only 30% yield — what could cause this?",
    ]
    for ex in examples:
        print(f"  {C.DIM}•{C.R}  {C.DIM}{ex}{C.R}")
    print()

    from chem_coworker import REGISTRY

    ui = TerminalUI()
    try:
        coworker = _init_coworker(model, provider, verbose, plan_mode, ui)
        public_tool_count = len(REGISTRY.filtered_names(llm_exposed_only=True))
        print(f"  {C.OK}✓{C.R}  Agent ready  {C.META}{public_tool_count} public tools available{C.R}")
        if plan_mode:
            print(f"  {C.META}Plan approval ON{C.R}")
    except Exception as exc:
        print(f"  {C.ERR}✗{C.R}  Failed to initialize: {exc}")
        sys.exit(1)

    session = _create_repl_session(coworker, model, provider, verbose, plan_mode, ui, REGISTRY)
    command_registry = build_default_command_registry()
    session.command_registry = command_registry

    while True:
        print()
        try:
            query = input(f"  {C.PROMPT}>{C.R} ").strip()
        except (KeyboardInterrupt, EOFError):
            print(f"\n\n  {C.META}Goodbye!{C.R}")
            break

        if not query:
            continue

        command_outcome = command_registry.dispatch(query, session)
        if command_outcome == COMMAND_EXIT:
            break
        if command_outcome == COMMAND_HANDLED:
            continue

        ui.reset_for_query()
        try:
            query_for_agent = _apply_condition_mode_hint(query, session.condition_mode)
            response, session.history = session.coworker.chat(query_for_agent, session.history)
        except Exception as exc:
            ui.stop_transient_ui()
            print(f"\n  {C.ERR}✗{C.R}  {exc}")
            continue

        if response.compacted:
            print(f"  {C.META}↺ Conversation history compacted{C.R}")

        ui.render_response(response, verbose=session.verbose)


if __name__ == "__main__":
    main()
