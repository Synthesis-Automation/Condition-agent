"""REPL command registry for the terminal CLI."""

from __future__ import annotations

import shlex
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Sequence

from .config import save_config
from .models import select_model_interactive
from .session_store import list_sessions, load_session, save_session
from .style import C


COMMAND_HANDLED = "handled"
COMMAND_EXIT = "exit"
COMMAND_UNHANDLED = "unhandled"


@dataclass
class ReplSession:
    """Mutable REPL state shared by command handlers."""

    coworker: Any
    model: str
    provider: str
    verbose: bool
    plan_mode: bool
    condition_mode: str
    history: List[Dict]
    tool_registry: Any
    ui: Any
    create_coworker: Callable[[str, str, bool, bool], Any]
    save_default_config: Callable[[str, str], None] | None = None
    command_registry: "CommandRegistry | None" = None


CommandHandler = Callable[[ReplSession, "CommandRegistry", str, List[str]], str]


@dataclass(frozen=True)
class Command:
    """One REPL command and its aliases."""

    name: str
    aliases: Sequence[str]
    help: str
    handler: CommandHandler


@dataclass
class CommandRegistry:
    """Exact-match command registry with slash-command parsing support."""

    _commands: List[Command] = field(default_factory=list)
    _aliases: Dict[str, Command] = field(default_factory=dict)

    def register(self, command: Command) -> "CommandRegistry":
        self._commands.append(command)
        for alias in command.aliases:
            self._aliases[alias.lower()] = command
        self._aliases[command.name.lower()] = command
        return self

    def dispatch(self, raw: str, session: ReplSession) -> str:
        text = raw.strip()
        if not text:
            return COMMAND_HANDLED

        cmd = self._aliases.get(text.lower())
        if cmd is not None:
            return self._invoke(cmd, session, text)

        if text.startswith("/"):
            try:
                parts = shlex.split(text)
            except ValueError as exc:
                print(f"  {C.ERR}✗{C.R}  {exc}")
                return COMMAND_HANDLED

            if not parts:
                return COMMAND_HANDLED
            cmd = self._aliases.get(parts[0].lower())
            if cmd is not None:
                return self._invoke(cmd, session, text, parts[1:])

        return COMMAND_UNHANDLED

    def _invoke(
        self,
        command: Command,
        session: ReplSession,
        raw: str,
        args: List[str] | None = None,
    ) -> str:
        try:
            parts = args if args is not None else shlex.split(raw)
        except ValueError:
            parts = [raw]
        return command.handler(session, self, raw, parts)

    def help_lines(self) -> List[str]:
        seen = set()
        lines: List[str] = []
        for cmd in self._commands:
            if cmd.name in seen:
                continue
            seen.add(cmd.name)
            lines.append(f"  {C.DIM}{cmd.name:<16}{C.R} {cmd.help}")
        return lines


def build_default_command_registry() -> CommandRegistry:
    """Create the default Codex-style REPL command set."""
    reg = CommandRegistry()
    reg.register(Command("/help", ("/help", "help", "?"), "Show CLI commands", _cmd_help))
    reg.register(Command("/exit", ("exit", "quit", "q", "/exit"), "Exit the session", _cmd_exit))
    reg.register(Command("/clear", ("clear", "/clear"), "Clear conversation history", _cmd_clear))
    reg.register(Command("/tools", ("tools", "/tools"), "List available tools", _cmd_tools))
    reg.register(Command("/model", ("/model", "change model", "model"), "Switch model/provider", _cmd_model))
    reg.register(Command("/plan", ("/plan", "toggle plan"), "Toggle plan approval mode", _cmd_plan))
    reg.register(Command("/verbose", ("/verbose", "toggle verbose"), "Toggle verbose tool output", _cmd_verbose))
    reg.register(Command("/condmode", ("/condmode", "condmode"), "Condition mode: auto|full", _cmd_condmode))
    reg.register(Command("/settings", ("/settings", "settings"), "Show current settings", _cmd_settings))
    reg.register(Command("/compact", ("/compact", "compact"), "Compact conversation history", _cmd_compact))
    reg.register(Command("/history", ("/history", "history"), "Show conversation history stats", _cmd_history))
    reg.register(Command("/session", ("/session",), "Session save/load/list/new", _cmd_session))
    return reg


def _cmd_help(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    print()
    print(f"  {C.LABEL}◆ Commands{C.R}")
    for line in registry.help_lines():
        print(line)
    print(f"  {C.DIM}Tip: /tools and /help support slash-command style like coding assistants.{C.R}")
    return COMMAND_HANDLED


def _cmd_exit(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    print(f"  {C.META}Goodbye!{C.R}")
    return COMMAND_EXIT


def _cmd_clear(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    session.history = []
    print(f"  {C.META}History cleared.{C.R}")
    return COMMAND_HANDLED


def _cmd_tools(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    sub = args[0].lower() if args else "public"
    if sub not in {"list", "public"}:
        print(f"  {C.WARN}⚠{C.R}  Unsupported /tools subcommand. Use `/tools` or `/tools public`.")
        return COMMAND_HANDLED
    print()
    print(session.tool_registry.describe_tools(llm_exposed_only=True))
    return COMMAND_HANDLED


def _cmd_model(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    selected = select_model_interactive()
    new_model = selected["name"]
    new_provider = selected["provider"]
    try:
        new_coworker = session.create_coworker(new_model, new_provider, session.verbose, session.plan_mode)
    except Exception as exc:
        print(f"  {C.ERR}✗{C.R}  Failed to switch model: {exc}")
        return COMMAND_HANDLED

    session.coworker = new_coworker
    session.model = new_model
    session.provider = new_provider
    if session.save_default_config is not None:
        session.save_default_config(new_model, new_provider)
    else:
        save_config(new_model, new_provider)
    provider_color = C.CYAN if new_provider == "openai" else C.MAGENTA
    print(f"  {C.OK}✓{C.R}  Switched to {C.BOLD}{new_model}{C.R}  {provider_color}{new_provider}{C.R}")
    print(f"  {C.META}Config saved. History preserved.{C.R}")
    return COMMAND_HANDLED


def _cmd_plan(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    session.plan_mode = not session.plan_mode
    session.coworker.set_plan_callback(session.ui.show_plan_and_confirm if session.plan_mode else None)
    state = f"{C.OK}ON{C.R}" if session.plan_mode else f"{C.DIM}off{C.R}"
    print(f"  {C.META}Plan approval{C.R}  {state}")
    return COMMAND_HANDLED


def _cmd_verbose(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    session.verbose = not session.verbose
    session.coworker.set_verbose(session.verbose)
    state = f"{C.OK}ON{C.R}" if session.verbose else f"{C.DIM}off{C.R}"
    print(f"  {C.META}Verbose{C.R}  {state}")
    return COMMAND_HANDLED


def _cmd_settings(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    session.ui.print_settings(
        session.model,
        session.provider,
        session.verbose,
        session.plan_mode,
        condition_mode=session.condition_mode,
    )
    return COMMAND_HANDLED


def _cmd_condmode(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    accepted = {"auto", "full", "balanced"}
    if not args:
        current = session.condition_mode
        print(f"  {C.META}Condition mode:{C.R} {C.BOLD}{current}{C.R}")
        print(f"  {C.DIM}Usage:{C.R} /condmode auto|full")
        return COMMAND_HANDLED

    mode = str(args[0]).strip().lower()
    if mode not in accepted:
        print(f"  {C.WARN}⚠{C.R}  Invalid mode: {mode}. Use auto or full.")
        return COMMAND_HANDLED
    if mode == "balanced":
        mode = "full"

    session.condition_mode = mode
    labels = {
        "auto": "No forced condition strategy (agent default behavior).",
        "full": "Bias condition requests toward full 4-mode condition analysis.",
    }
    print(f"  {C.OK}✓{C.R}  Condition mode set to {C.BOLD}{mode}{C.R}")
    print(f"  {C.META}{labels[mode]}{C.R}")
    return COMMAND_HANDLED


def _cmd_compact(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    if len(session.history) < 2:
        print(f"  {C.META}Nothing to compact (history is empty).{C.R}")
        return COMMAND_HANDLED
    session.history = session.coworker.compact_history(session.history)
    print(f"  {C.META}↺ History compacted to {len(session.history)} messages.{C.R}")
    return COMMAND_HANDLED


def _cmd_history(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    turns = len(session.history) // 2
    print(f"  {C.META}History:{C.R} {len(session.history)} messages ({turns} turns)")
    return COMMAND_HANDLED


def _cmd_session(session: ReplSession, registry: CommandRegistry, raw: str, args: List[str]) -> str:  # noqa: ARG001
    if not args:
        print(f"  {C.DIM}Usage:{C.R} /session save [name] | load <name> | list | new")
        return COMMAND_HANDLED

    sub = args[0].lower()

    if sub == "save":
        name = args[1] if len(args) >= 2 else None
        try:
            path = save_session(
                name=name,
                history=session.history,
                model=session.model,
                provider=session.provider,
                verbose=session.verbose,
                plan_mode=session.plan_mode,
                condition_mode=session.condition_mode,
            )
        except Exception as exc:
            print(f"  {C.ERR}✗{C.R}  Failed to save session: {exc}")
            return COMMAND_HANDLED
        print(f"  {C.OK}✓{C.R}  Session saved: {C.DIM}{path}{C.R}")
        return COMMAND_HANDLED

    if sub == "list":
        files = list_sessions()
        if not files:
            print(f"  {C.META}No saved sessions found.{C.R}")
            return COMMAND_HANDLED
        print()
        print(f"  {C.LABEL}◆ Saved Sessions{C.R}")
        for path in files:
            print(f"  {C.DIM}{path.stem}{C.R}  {C.META}{path}{C.R}")
        return COMMAND_HANDLED

    if sub == "new":
        session.history = []
        print(f"  {C.META}Started a new session (history cleared).{C.R}")
        return COMMAND_HANDLED

    if sub == "load":
        if len(args) < 2:
            print(f"  {C.DIM}Usage:{C.R} /session load <name>")
            return COMMAND_HANDLED
        name = args[1]
        try:
            data = load_session(name)
        except Exception as exc:
            print(f"  {C.ERR}✗{C.R}  Failed to load session: {exc}")
            return COMMAND_HANDLED

        target_model = str(data.get("model") or session.model)
        target_provider = str(data.get("provider") or session.provider)
        if (target_model, target_provider) != (session.model, session.provider):
            try:
                session.coworker = session.create_coworker(
                    target_model, target_provider, session.verbose, session.plan_mode
                )
                session.model = target_model
                session.provider = target_provider
                if session.save_default_config is not None:
                    session.save_default_config(target_model, target_provider)
            except Exception as exc:
                print(
                    f"  {C.WARN}⚠{C.R}  Loaded history, but could not switch to saved model "
                    f"{target_model}/{target_provider}: {exc}"
                )

        session.history = list(data.get("history", []))
        session.condition_mode = str(data.get("condition_mode") or session.condition_mode or "auto").lower()
        if session.condition_mode == "balanced":
            session.condition_mode = "full"
        if session.condition_mode not in {"auto", "full"}:
            session.condition_mode = "auto"
        print(
            f"  {C.OK}✓{C.R}  Session loaded: {C.DIM}{name}{C.R}  "
            f"{C.META}{len(session.history)} messages{C.R}"
        )
        print(f"  {C.META}Condition mode:{C.R} {C.BOLD}{session.condition_mode}{C.R}")
        return COMMAND_HANDLED

    print(f"  {C.WARN}⚠{C.R}  Unknown /session subcommand: {sub}")
    print(f"  {C.DIM}Usage:{C.R} /session save [name] | load <name> | list | new")
    return COMMAND_HANDLED
