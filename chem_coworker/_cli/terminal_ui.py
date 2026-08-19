"""Rich, prompt-toolkit frontend for the interactive condition coworker."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, Optional

from prompt_toolkit import PromptSession
from prompt_toolkit.auto_suggest import AutoSuggestFromHistory
from prompt_toolkit.completion import Completer, Completion
from prompt_toolkit.document import Document
from prompt_toolkit.formatted_text import HTML
from prompt_toolkit.history import FileHistory, InMemoryHistory
from prompt_toolkit.key_binding import KeyBindings
from prompt_toolkit.styles import Style
from rich import box
from rich.console import Console
from rich.json import JSON
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from condition_recommender import available_ranking_profiles

from chem_coworker.contracts import ConditionResponse

from .interactive import InteractiveSession, InteractiveSettings, _InteractiveCoworker


HISTORY_PATH = Path.home() / ".chemcoworker" / "history"

COMMANDS = {
    "/help": "show commands and keyboard shortcuts",
    "/settings": "show recommendation settings",
    "/top-k": "set the recommendation count",
    "/minimum": "set the retrieval pool minimum",
    "/profile": "list or select a ranking profile",
    "/json": "toggle full JSON output",
    "/last": "show the previous result",
    "/save": "save the previous result",
    "/clear": "clear the terminal",
    "/quit": "exit the coworker",
}


def can_use_terminal_ui() -> bool:
    """Return whether both standard streams are interactive terminals."""

    return bool(sys.stdin.isatty() and sys.stdout.isatty())


class SlashCommandCompleter(Completer):
    """Complete slash commands and their finite argument vocabularies."""

    def __init__(self, profile_ids: Iterable[str]) -> None:
        self.profile_ids = tuple(profile_ids)

    def get_completions(self, document: Document, complete_event):  # type: ignore[no-untyped-def]
        text = document.text_before_cursor
        if not text.startswith("/"):
            return
        command, separator, argument = text.partition(" ")
        if not separator:
            for candidate, description in COMMANDS.items():
                if candidate.startswith(command.casefold()):
                    yield Completion(
                        candidate,
                        start_position=-len(command),
                        display_meta=description,
                    )
            return
        options: tuple[str, ...] = ()
        if command.casefold() == "/profile":
            options = self.profile_ids
        elif command.casefold() == "/json":
            options = ("on", "off")
        elif command.casefold() == "/minimum":
            options = ("auto",)
        for candidate in options:
            if candidate.startswith(argument):
                yield Completion(candidate, start_position=-len(argument))


class RichResponseRenderer:
    """Render typed recommendation results without changing their content."""

    _RECIPE_GROUPS = (
        ("catalysts", "catalyst"),
        ("ligands", "ligand"),
        ("bases", "base"),
        ("condensation_agents", "coupling reagent"),
        ("oxidants", "oxidant"),
        ("reductants", "reductant"),
        ("acids", "acid"),
        ("additives", "additive"),
        ("solvents", "solvent"),
        ("other_components", "other"),
    )

    def __init__(self, console: Console) -> None:
        self.console = console

    def __call__(self, response: ConditionResponse, as_json: bool) -> None:
        if as_json:
            self.console.print(JSON.from_data(response.to_dict()))
            return
        result = response.result
        if not result.valid:
            self.console.print(
                Panel(
                    Text(result.error or "Invalid reaction", style="bold red"),
                    title="Recommendation unavailable",
                    border_style="red",
                )
            )
            self._warnings(result.warnings)
            return

        label = result.reaction_label or {}
        summary = Table.grid(padding=(0, 2))
        summary.add_column(style="bold cyan")
        summary.add_column()
        summary.add_row(
            "Reaction",
            Text(
                str(label.get("concise") or label.get("detailed") or "Unlabeled")
            ),
        )
        summary.add_row(
            "Transformation", Text(result.transformation_class or "unresolved")
        )
        summary.add_row("Named family", Text(result.named_family or "unassigned"))
        summary.add_row("Retrieval", Text(result.retrieval_level or "none"))
        summary.add_row("Candidates", str(result.candidate_count))
        self.console.print(
            Panel(summary, title="Reaction analysis", border_style="cyan")
        )

        if not result.recommendations:
            self.console.print(
                Panel(
                    "No compatible recipe met the recommendation criteria.",
                    border_style="yellow",
                )
            )
            self._warnings(result.warnings)
            return

        table = Table(
            title="Recommended conditions",
            box=box.ROUNDED,
            header_style="bold magenta",
            show_lines=True,
        )
        table.add_column("#", justify="right", style="bold")
        table.add_column("Resolved recipe", ratio=4)
        table.add_column("Score", justify="right")
        table.add_column("Yield", justify="right")
        table.add_column("Refs", justify="right")
        for item in result.recommendations:
            recipe = self._recipe(item.resolved_recipe, item.recipe_core_id)
            expected_yield = (
                f"{item.expected_yield_pct:.1f}%"
                if item.expected_yield_pct is not None
                else "—"
            )
            table.add_row(
                str(item.rank),
                Text(recipe),
                f"{item.score:.3f}",
                expected_yield,
                str(item.reference_support),
            )
        self.console.print(table)

        for item in result.recommendations:
            details = Text()
            if item.explanation:
                details.append("Evidence\n", style="bold green")
                details.append("\n".join(f"• {value}" for value in item.explanation))
            if item.cautions:
                if details:
                    details.append("\n\n")
                details.append("Cautions\n", style="bold yellow")
                details.append("\n".join(f"• {value}" for value in item.cautions))
            if item.precedent_reaction_ids:
                if details:
                    details.append("\n\n")
                details.append("Precedents\n", style="bold cyan")
                details.append(", ".join(item.precedent_reaction_ids))
            if details:
                self.console.print(
                    Panel(
                        details,
                        title=f"Recipe {item.rank} details",
                        border_style="dim",
                    )
                )
        self._warnings(result.warnings)

    def _warnings(self, warnings: tuple[str, ...]) -> None:
        if warnings:
            body = Text()
            body.append("\n".join(f"• {warning}" for warning in warnings))
            self.console.print(
                Panel(
                    body,
                    title="Warnings",
                    border_style="yellow",
                )
            )

    def _recipe(self, recipe: dict[str, object], fallback: str) -> str:
        values = []
        for key, role in self._RECIPE_GROUPS:
            components = recipe.get(key) or ()
            if not isinstance(components, (list, tuple)):
                continue
            for component in components:
                if not isinstance(component, dict):
                    continue
                name = (
                    component.get("canonical_name")
                    or component.get("raw_identifier")
                    or component.get("substance_id")
                )
                if name:
                    values.append(f"{name} ({role})")
        return "\n".join(values) if values else fallback


def _key_bindings() -> KeyBindings:
    bindings = KeyBindings()

    @bindings.add("enter")
    def submit(event) -> None:  # type: ignore[no-untyped-def]
        event.current_buffer.validate_and_handle()

    @bindings.add("escape", "enter")
    def newline(event) -> None:  # type: ignore[no-untyped-def]
        event.current_buffer.insert_text("\n")

    return bindings


def run_terminal_ui(
    coworker: _InteractiveCoworker,
    *,
    settings: Optional[InteractiveSettings] = None,
    initial_reaction: Optional[str] = None,
    persistent_history: bool = True,
) -> int:
    """Run the rich interactive terminal application."""

    console = Console(highlight=False)
    active_settings = settings or InteractiveSettings()
    profile_ids = tuple(
        str(item["profile_id"]) for item in available_ranking_profiles()
    )
    if persistent_history:
        try:
            HISTORY_PATH.parent.mkdir(parents=True, exist_ok=True)
            history = FileHistory(str(HISTORY_PATH))
        except OSError as exc:
            message = Text(
                f"History is unavailable ({exc}); using session-only history.",
                style="yellow",
            )
            console.print(message)
            history = InMemoryHistory()
    else:
        history = InMemoryHistory()
    prompt = PromptSession(
        history=history,
        auto_suggest=AutoSuggestFromHistory(),
        completer=SlashCommandCompleter(profile_ids),
        complete_while_typing=False,
        key_bindings=_key_bindings(),
        multiline=True,
        style=Style.from_dict(
            {
                "prompt": "ansicyan bold",
                "bottom-toolbar": "bg:#222222 #aaaaaa",
                "completion-menu.completion.current": "bg:#005f87 #ffffff",
            }
        ),
    )

    def read_input(message: str) -> str:
        prompt_label = "source › " if message.startswith("source") else (
            "custom › " if message.startswith("custom") else "❯ "
        )
        return prompt.prompt(
            HTML(f"<prompt>{prompt_label}</prompt>"),
            bottom_toolbar=lambda: HTML(
                " <b>Enter</b> send  <b>Alt+Enter</b> newline  "
                f"<b>profile</b> {active_settings.ranking_profile}  "
                f"<b>top-k</b> {active_settings.top_k}  <b>/help</b> commands "
            ),
        )

    def write_output(value: str) -> None:
        console.print(Text(value))

    session = InteractiveSession(
        coworker,
        settings=active_settings,
        input_fn=read_input,
        output_fn=write_output,
        response_renderer=RichResponseRenderer(console),
        status_fn=lambda message: console.status(message, spinner="dots"),
        clear_fn=console.clear,
    )
    console.print(
        Panel.fit(
            "[bold cyan]Condition Coworker[/]\n"
            "Structure-first recommendation with auditable evidence\n\n"
            "Type a reaction SMILES or [bold]/help[/]",
            border_style="cyan",
        )
    )
    return session.run(initial_reaction, show_welcome=False)


__all__ = [
    "HISTORY_PATH",
    "RichResponseRenderer",
    "SlashCommandCompleter",
    "can_use_terminal_ui",
    "run_terminal_ui",
]
