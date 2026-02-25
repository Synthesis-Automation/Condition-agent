"""Terminal UI rendering and event-bus observer wiring."""

from __future__ import annotations

import dataclasses
import json
import sys
import threading
import time
from typing import TYPE_CHECKING, Any, Dict, Optional

from chem_coworker.event_bus import ChemEvent, EventBus

from .config import CONFIG_PATH
from .style import C, SEP_FAT, label

if TYPE_CHECKING:
    from chem_coworker.response import ChemResponse


class Spinner:
    """Braille spinner for phase/progress display."""

    _FRAMES = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]

    def __init__(self, lock: threading.Lock, message: str = "Thinking") -> None:
        self._lock = lock
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
        with self._lock:
            sys.stdout.write("\r\033[2K")
            sys.stdout.flush()

    def _spin(self) -> None:
        idx = 0
        while not self._stop.is_set():
            frame = self._FRAMES[idx % len(self._FRAMES)]
            with self._lock:
                sys.stdout.write(
                    f"\r\033[2K  {C.YELLOW}{frame}{C.R} {C.DIM}{self._message}…{C.R}"
                )
                sys.stdout.flush()
            time.sleep(0.1)
            idx += 1


class _LineWriter:
    """Line-buffered writer for streamed answer output."""

    def __init__(self) -> None:
        self._buf = ""

    def write(self, token: str) -> None:
        self._buf += token
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            sys.stdout.write(f"  {line}\n")
        sys.stdout.flush()

    def flush_remaining(self) -> None:
        if self._buf:
            sys.stdout.write(f"  {self._buf}")
            sys.stdout.flush()
            self._buf = ""

    def reset(self) -> None:
        self._buf = ""


class TerminalUI:
    """Session-scoped terminal UI state + renderers + event observers."""

    _PHASE_LABELS = {
        "reason_start": "Reasoning",
        "observe_start": "Observing",
        "synth_start": "Synthesizing",
        "compact_start": "Compacting history",
    }

    def __init__(self) -> None:
        self.stdout_lock = threading.Lock()
        self.running_tools: set[str] = set()
        self.phase_spinner: Optional[Spinner] = None
        self.stream_state: Dict[str, Any] = {
            "first_chunk": True,
            "pre_synth_info": {},
            "writer": _LineWriter(),
        }

    def reset_for_query(self) -> None:
        self.running_tools.clear()
        self.stream_state["first_chunk"] = True
        self.stream_state["pre_synth_info"] = {}
        self.stream_state["writer"].reset()

    def stop_transient_ui(self) -> None:
        if self.phase_spinner is not None:
            self.phase_spinner.stop()
            self.phase_spinner = None

    def spinner(self, message: str) -> Spinner:
        return Spinner(self.stdout_lock, message)

    def build_event_bus(self) -> EventBus:
        """Wire this UI instance to a fresh EventBus."""
        bus = EventBus()
        bus.subscribe(ChemEvent.TOOL_START, lambda tool_name, **_: self._progress("start", tool_name, 0.0))
        bus.subscribe(ChemEvent.TOOL_DONE, lambda tool_name, elapsed_s, **_: self._progress("done", tool_name, elapsed_s))
        bus.subscribe(ChemEvent.TOOL_ERROR, lambda tool_name, elapsed_s, **_: self._progress("error", tool_name, elapsed_s))

        bus.subscribe(ChemEvent.PHASE_START, lambda phase, **_: self._phase(f"{phase}_start"))
        bus.subscribe(ChemEvent.PHASE_END, lambda phase, **_: self._phase(f"{phase}_done"))
        bus.subscribe(ChemEvent.COMPACT_START, lambda **_: self._phase("compact_start"))
        bus.subscribe(ChemEvent.COMPACT_END, lambda **_: self._phase("compact_done"))

        bus.subscribe(
            ChemEvent.PRE_SYNTH,
            lambda hypothesis, confidence, rationale, tools_called, **_:
                self._pre_synth_cb(hypothesis, confidence, rationale, tools_called),
        )
        bus.subscribe(ChemEvent.STREAM_TOKEN, lambda token, **_: self._stream_token(token))
        return bus

    def _progress(self, event: str, name: str, elapsed: float) -> None:
        if event == "start":
            self.running_tools.add(name)
            with self.stdout_lock:
                sys.stdout.write(
                    f"\r\033[2K  {C.YELLOW}⠋{C.R} {C.DIM}{' · '.join(sorted(self.running_tools))}…{C.R}"
                )
                sys.stdout.flush()
            return

        self.running_tools.discard(name)
        icon = f"{C.OK}✓" if event == "done" else f"{C.ERR}✗"
        with self.stdout_lock:
            sys.stdout.write(
                f"\r\033[2K  {icon}{C.R}  {C.TOOL}{name}{C.R}  {C.META}{elapsed:.1f}s{C.R}\n"
            )
            sys.stdout.flush()

    def _phase(self, phase: str) -> None:
        if phase.endswith("_start") and phase in self._PHASE_LABELS:
            self.phase_spinner = Spinner(self.stdout_lock, self._PHASE_LABELS[phase])
            self.phase_spinner.start()
        elif phase.endswith("_done") and self.phase_spinner is not None:
            self.phase_spinner.stop()
            self.phase_spinner = None

    def _pre_synth_cb(
        self,
        hypothesis: str,
        confidence: float,
        rationale: str,
        tools_called: list,
    ) -> None:
        self.stream_state["pre_synth_info"] = {
            "hypothesis": hypothesis,
            "confidence": confidence,
            "rationale": rationale,
            "tools_called": tools_called,
        }

    def _stream_token(self, token: str) -> None:
        if self.stream_state["first_chunk"]:
            self.stream_state["first_chunk"] = False
            if self.phase_spinner is not None:
                self.phase_spinner.stop()
                self.phase_spinner = None
            self._print_pre_answer_chrome(self.stream_state["pre_synth_info"])
            print(SEP_FAT)
            self.stream_state["writer"].reset()

        self.stream_state["writer"].write(token)

    def _print_pre_answer_chrome(self, ctx: dict) -> None:
        hypothesis = ctx.get("hypothesis", "")
        confidence = ctx.get("confidence", 0.0)
        rationale = ctx.get("rationale", "")
        tools_called = ctx.get("tools_called", [])

        print()
        if hypothesis or rationale:
            conf_str = f"{C.CONF}[{confidence:.0%}]{C.R}" if confidence else ""
            print(label("◆", "Hypothesis") + f"  {conf_str}")
            if hypothesis:
                print(f"  {C.HYPO}{hypothesis}{C.R}")
            if rationale:
                print(f"  {C.META}{rationale}{C.R}")
            print()

        if tools_called:
            arrow = f"  {C.META}→{C.R}  "
            tool_chain = arrow.join(f"{C.TOOL}{t}{C.R}" for t in tools_called)
            print(label("⎿", "Tools"))
            print(f"  {tool_chain}")
        else:
            print(f"  {C.META}⎿ (no tools called — answered from LLM knowledge){C.R}")
        print()

    def render_response(self, response: "ChemResponse", verbose: bool = False) -> None:
        """Render the final response (including streamed or non-streamed modes)."""
        if response.streamed:
            self.stream_state["writer"].flush_remaining()
            print()
            print(SEP_FAT)
            if verbose and response.tool_results:
                print()
                self._print_tool_details(response.tool_results)
        else:
            print()
            if response.hypothesis or response.plan_rationale:
                conf_str = f"{C.CONF}[{response.confidence:.0%}]{C.R}" if response.confidence else ""
                print(label("◆", "Hypothesis") + f"  {conf_str}")
                if response.hypothesis:
                    print(f"  {C.HYPO}{response.hypothesis}{C.R}")
                if response.plan_rationale:
                    print(f"  {C.META}{response.plan_rationale}{C.R}")
                print()

            if response.tools_called:
                arrow = f"  {C.META}→{C.R}  "
                tool_chain = arrow.join(f"{C.TOOL}{t}{C.R}" for t in response.tools_called)
                print(label("⎿", "Tools"))
                print(f"  {tool_chain}")
                if verbose and response.tool_results:
                    print()
                    self._print_tool_details(response.tool_results)
            else:
                print(f"  {C.META}⎿ (no tools called — answered from LLM knowledge){C.R}")

            print()
            print(SEP_FAT)
            for line in response.answer.strip().splitlines():
                print(f"  {line}" if line.strip() else "")
            print(SEP_FAT)

        parts = [
            response.model,
            f"{response.elapsed_s:.1f}s",
            f"{response.llm_calls} LLM calls",
            f"{len(response.tools_called)} tools",
        ]
        print(f"  {C.META}{' · '.join(parts)}{C.R}")

        if response.warnings:
            print()
            for w in response.warnings:
                print(f"  {C.WARN}⚠  {w}{C.R}")

    def _print_tool_details(self, tool_results: Dict[str, Any]) -> None:
        for name, result in tool_results.items():
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

    def show_plan_and_confirm(self, plan):
        """Display a plan and pause for approval."""
        from chem_coworker.plan import PlanRejected

        print()
        conf_str = f"{C.CONF}[{plan.confidence:.0%}]{C.R}" if plan.confidence else ""
        print(label("◆", "Proposed Plan") + f"  {conf_str}")
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
            return self._drop_tool_from_plan(plan, raw[5:].strip())
        print()
        return plan

    def _drop_tool_from_plan(self, plan, name: str):
        new_groups = [[c for c in group if c.name != name] for group in plan.groups]
        new_groups = [g for g in new_groups if g]
        return dataclasses.replace(plan, groups=new_groups)

    def print_settings(self, model: str, provider: str, verbose: bool, plan_mode: bool) -> None:
        provider_color = C.CYAN if provider == "openai" else C.MAGENTA
        print()
        print(label("◆", "Current Settings"))
        print(f"  {C.META}Model   {C.R}  {C.BOLD}{model}{C.R}  {provider_color}{provider}{C.R}")
        print(f"  {C.META}Verbose {C.R}  {C.OK if verbose else C.DIM}{'on' if verbose else 'off'}{C.R}")
        print(f"  {C.META}Plan    {C.R}  {C.OK if plan_mode else C.DIM}{'on' if plan_mode else 'off'}{C.R}")
        print(f"  {C.META}Config  {C.R}  {C.DIM}{CONFIG_PATH}{C.R}")
        print()
        print(f"  {C.DIM}Commands: /help · /model · /plan · /verbose · /tools · /compact · /session · /settings{C.R}")
