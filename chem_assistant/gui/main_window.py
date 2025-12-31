"""Main Qt window for the ChemTools assistant."""

from __future__ import annotations

import html
import json
import re
import shlex
import tempfile
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlparse
from uuid import uuid4

from PyQt6.QtCore import QObject, QRunnable, Qt, QThreadPool, pyqtSignal, QTimer
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QStatusBar,
    QTextEdit,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)
from PyQt6.QtGui import QFont, QPixmap

from langchain_core.messages import BaseMessage

from chem_assistant.chemtools_agent import ChemToolsAgent
from chem_assistant.constraint_parser import (
    ConstraintSpec,
    build_constraint_spec,
    format_constraints_for_prompt,
    merge_specs,
)
from chem_assistant.gui.dialogs import (
    ConstraintDialog,
    LLMConfigDialog,
)
from chemtools.visualization import render_molecule_image, render_reaction_image
from chemtools.util.rdkit_helpers import parse_smiles, rdkit_available
from chemtools.featurizers.reaction_detection import detect_reaction_types

IMAGE_MARKUP = re.compile(
    r"\[\[(reaction|molecule)_image:(.+?)\]\]", re.IGNORECASE | re.DOTALL
)
MARKDOWN_IMAGE = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")
IMAGE_COMMAND_PATTERNS = [
    r"^(?:/)?image\s+(reaction|molecule|compound)\s*[:=]\s*(.+)$",
    r"^(?:/)?image\s+(?:for\s+)?(reaction|molecule|compound)\s*[:=]\s*(.+)$",
    r"^(?:show|display)\s+(?:the\s+)?image(?:\s+of|\s+for)?\s*(?:a|the)?\s*(reaction|molecule|compound)?[:=]?\s*(.+)$",
    r"^(?:show|display)\s+(?:me\s+)?(?:an?\s+|the\s+)?image\s*(?:of|for)?\s*(reaction|molecule|compound)?[:=]?\s*(.+)$",
]
URL_SMILES_KEYS = ("smiles", "model", "structure", "mol", "target")
REACTION_SMILES = re.compile(
    r"([A-Za-z0-9@+\-\[\]\(\)=#$\.]+>>[A-Za-z0-9@+\-\[\]\(\)=#$\.]+)"
)
CONDITION_KEYWORDS = (
    "condition",
    "conditions",
    "recommend",
    "recommendation",
    "setup",
    "protocol",
    "hte",
    "catalyst",
    "ligand",
    "base",
    "solvent",
)


class WorkerSignals(QObject):
    """Signals for background tasks."""

    finished = pyqtSignal()
    error = pyqtSignal(str)
    result = pyqtSignal(object)
    progress = pyqtSignal(str)


class Worker(QRunnable):
    """Run blocking work on a background thread."""

    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    def run(self) -> None:
        try:
            # Inject progress callback if the function accepts it
            if "progress_callback" in self.fn.__code__.co_varnames:
                self.kwargs["progress_callback"] = self.signals.progress.emit
            result = self.fn(*self.args, **self.kwargs)
        except Exception as exc:  # pragma: no cover - UI feedback
            trace = "".join(traceback.format_exception(exc))
            self.signals.error.emit(trace)
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()


class ChatInput(QTextEdit):
    """Text edit that emits a signal on Ctrl+Enter."""

    sendRequested = pyqtSignal()

    def keyPressEvent(self, event):  # type: ignore[override]
        if (
            event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter)
            and event.modifiers() & Qt.KeyboardModifier.ControlModifier
        ):
            event.accept()
            self.sendRequested.emit()
            return
        super().keyPressEvent(event)


class ImagePreviewWindow(QDialog):
    """Standalone window that renders a chemistry image."""

    def __init__(
        self,
        pixmap: QPixmap,
        caption: str,
        *,
        parent: Optional[QWidget] = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(caption or "ChemTools Image Preview")
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setWindowModality(Qt.WindowModality.NonModal)

        self._original_pixmap = pixmap

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)

        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.image_label.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        layout.addWidget(self.image_label)

        if caption:
            caption_label = QLabel(caption)
            caption_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            caption_label.setWordWrap(True)
            caption_label.setStyleSheet("font-weight: 600; padding-top: 6px;")
            layout.addWidget(caption_label)

        self.resize(
            max(420, pixmap.width() + 64),
            max(360, pixmap.height() + 96),
        )
        self._update_pixmap()

    def resizeEvent(self, event):  # type: ignore[override]
        super().resizeEvent(event)
        self._update_pixmap()

    def _update_pixmap(self) -> None:
        if self._original_pixmap.isNull():
            return
        target_size = self.image_label.size()
        if not target_size.isValid():
            return
        scaled = self._original_pixmap.scaled(
            target_size,
            Qt.AspectRatioMode.KeepAspectRatio,
            Qt.TransformationMode.SmoothTransformation,
        )
        self.image_label.setPixmap(scaled)


class JsonViewerDialog(QDialog):
    """Simple dialog to view JSON content."""

    def __init__(self, data: Dict[str, Any], title: str = "JSON Viewer", parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(800, 600)
        layout = QVBoxLayout(self)
        
        self.text_edit = QTextEdit()
        self.text_edit.setReadOnly(True)
        self.text_edit.setPlainText(json.dumps(data, indent=2))
        
        # Use a monospace font
        font = QFont("Courier New")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(10)
        self.text_edit.setFont(font)
        
        layout.addWidget(self.text_edit)
        
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)


class ChemAssistantWindow(QMainWindow):
    """Main application window."""

    def __init__(self, *, startup_message: str | None = None) -> None:
        super().__init__()
        self.setWindowTitle("ChemTools Assistant")
        self.resize(1200, 800)
        self._apply_default_font()

        self.agent: Optional[ChemToolsAgent] = None
        self.llm_provider: Optional[str] = None
        self.llm_model: Optional[str] = None
        self.llm_temperature: float = 0.0
        self.proactive_enabled: bool = False
        self.proactive_top_k: int = 5
        self.proactive_max_protocols: int = 3
        self.proactive_build_protocols: bool = True
        self.history: List[BaseMessage] = []
        self.constraint_spec = ConstraintSpec()
        self.constraint_text = ""
        self.current_image_path: Optional[Path] = None
        self.last_details_path: Optional[str] = None
        self._image_windows: List[ImagePreviewWindow] = []
        self._startup_message = (startup_message or "").strip()

        self.thread_pool = QThreadPool()
        self._build_ui()
        self._update_constraint_summary()
        self.setAcceptDrops(True)
        self.statusBar().showMessage(self._startup_message or "Ready")
        if self._startup_message:
            self.append_log(self._startup_message)

    # ------------------------------------------------------------------ #
    # Styling helpers
    # ------------------------------------------------------------------ #

    def _apply_default_font(self) -> None:
        font = QFont(self.font())
        if font.pointSize() > 0:
            font.setPointSize(font.pointSize() + 2)
        else:
            font.setPointSize(12)
        self.setFont(font)

    # ------------------------------------------------------------------ #
    # UI construction
    # ------------------------------------------------------------------ #

    def _build_ui(self) -> None:
        central = QWidget()
        central_layout = QVBoxLayout(central)
        central_layout.setContentsMargins(8, 8, 8, 8)
        central_layout.setSpacing(8)

        self.chat_view = QTextEdit()
        self.chat_view.setReadOnly(True)
        self.chat_view.setAcceptRichText(True)
        self.chat_view.setPlaceholderText("Conversation history will appear here.")
        self.chat_view.setMinimumHeight(480)
        central_layout.addWidget(self.chat_view)

        # Input area
        self.input_edit = ChatInput()
        self.input_edit.setPlaceholderText(
            "Type a question, use /conditions <reaction_smiles>, or paste a reaction SMILES..."
        )
        self.input_edit.setMinimumHeight(60)
        self.input_edit.setMaximumHeight(140)
        self.input_edit.setStyleSheet(
            "QTextEdit { background-color: #30343b; color: #f5f6fa; border: 1px solid #555; }"
        )
        input_font = QFont(self.font())
        input_font.setPointSize(input_font.pointSize() + 1)
        self.input_edit.setFont(input_font)
        self.input_edit.sendRequested.connect(self.send_message)
        self.input_edit.textChanged.connect(self._on_input_changed)
        central_layout.addWidget(self.input_edit)

        self.validation_label = QLabel("")
        self.validation_label.setStyleSheet("color: #888; font-size: 10pt;")
        central_layout.addWidget(self.validation_label)

        button_row = QHBoxLayout()
        self.send_btn = QPushButton("Send")
        self.send_btn.clicked.connect(self.send_message)
        button_row.addWidget(self.send_btn)

        clear_btn = QPushButton("Clear Chat")
        clear_btn.clicked.connect(self.clear_chat)
        button_row.addWidget(clear_btn)

        constraint_btn = QPushButton("Manage Constraints")
        constraint_btn.clicked.connect(self.manage_constraints)
        button_row.addWidget(constraint_btn)

        llm_btn = QPushButton("LLM Mode")
        llm_btn.clicked.connect(self.configure_llm_mode)
        button_row.addWidget(llm_btn)

        self.view_results_btn = QPushButton("View Last Results")
        self.view_results_btn.clicked.connect(self.view_last_results)
        self.view_results_btn.setEnabled(False)
        button_row.addWidget(self.view_results_btn)

        button_row.addStretch()
        central_layout.addLayout(button_row)

        self.setCentralWidget(central)
        self.setStatusBar(QStatusBar())

    # ------------------------------------------------------------------ #
    # Chat handling
    # ------------------------------------------------------------------ #

    def append_chat(self, speaker: str, text: str) -> None:
        sanitized = html.escape(text.strip()).replace("\n", "<br>")
        if speaker.lower().startswith("you"):
            speaker_color = "#1f6feb"
            text_color = "#1f6feb"
        else:
            speaker_color = "#f0f4ff"
            text_color = "#f8f9fb"
        html_block = (
            f'<div style="margin-bottom:6px;">'
            f'<span style="color:{speaker_color}; font-weight:bold;">{html.escape(speaker)}:</span> '
            f'<span style="color:{text_color};">{sanitized}</span>'
            f"</div>"
        )
        self.chat_view.append(html_block)

    def append_log(self, text: str) -> None:
        message = text.strip()
        if not message:
            return
        self.statusBar().showMessage(message, 5000)

    def _on_input_changed(self) -> None:
        content = self.input_edit.toPlainText().strip()
        if not content:
            self.validation_label.setText("")
            return

        # Check for reaction SMILES
        reaction_smiles = self._extract_reaction_smiles(content)
        if reaction_smiles:
            if ">>" in reaction_smiles:
                parts = reaction_smiles.split(">>")
                if len(parts) == 2:
                    reactants, products = parts
                    r_mol = parse_smiles(reactants)
                    p_mol = parse_smiles(products)
                    if r_mol and p_mol:
                        # Try to detect reaction type
                        try:
                            detection = detect_reaction_types(reaction_smiles)
                            if detection.matches:
                                best = detection.matches[0]
                                self.validation_label.setText(f"✓ Valid Reaction: {best.name} detected")
                            else:
                                self.validation_label.setText("✓ Valid Reaction SMILES detected")
                        except Exception:
                            self.validation_label.setText("✓ Valid Reaction SMILES detected")
                        
                        self.validation_label.setStyleSheet("color: #4cd137;")
                        return
            elif parse_smiles(reaction_smiles):
                self.validation_label.setText("✓ Valid Molecule SMILES detected")
                self.validation_label.setStyleSheet("color: #4cd137;")
                return

        # If it looks like a command, don't show validation error
        if content.startswith("/"):
            self.validation_label.setText("")
            return

        # If it's long and doesn't look like SMILES, clear
        if len(content) > 10 and not any(c in content for c in "()[]=#123456"):
            self.validation_label.setText("")
            return

        self.validation_label.setText("")

    def send_message(self) -> None:
        content = self.input_edit.toPlainText().strip()
        if not content:
            return
        self.input_edit.clear()
        if self.handle_local_command(content):
            return
        if self._maybe_start_conditions_workflow(content):
            return
        self.append_chat("You", content)
        if not self._ensure_agent():
            self.append_log("Agent unavailable; use local commands.")
            return
        self.statusBar().showMessage("Agent is thinking... |")
        self.set_input_enabled(False)
        self.start_spinner("Agent is thinking")

        worker = Worker(
            self.agent.chat,  # type: ignore[union-attr]
            content,
            self.history,
            15,
            constraints=self.constraint_spec,
        )
        worker.signals.result.connect(self.handle_agent_result)
        worker.signals.error.connect(self.handle_agent_error)
        worker.signals.finished.connect(self.reset_after_task)
        self.thread_pool.start(worker)

    def handle_local_command(self, content: str) -> bool:
        command_text = content.strip()
        if not command_text:
            return False

        parsed_tokens = None
        if command_text.startswith("/"):
            try:
                parsed_tokens = shlex.split(command_text)
            except ValueError:
                parsed_tokens = [command_text]
        if parsed_tokens:
            cmd = parsed_tokens[0].lstrip("/").lower()
            args = parsed_tokens[1:]
            if cmd in {"taxonomy", "tax"}:
                self._handle_taxonomy_command()
                return True
            if cmd in {"terms", "term"}:
                self._handle_terms_command(cmd, args)
                return True
            if cmd in {"conditions", "condition", "cond", "auto", "recommend"}:
                self._handle_conditions_command(command_text, args)
                return True

        parsed = self._extract_image_request(command_text)
        if parsed:
            target, smiles = parsed
            self._execute_image_request(target, smiles, log_to_chat=True)
            return True
        return False

    def _extract_reaction_smiles(self, content: str) -> Optional[str]:
        """Return a best-effort reaction SMILES substring (reactants>>products)."""
        text = (content or "").strip()
        if not text or ">>" not in text:
            return None

        # Strip common wrappers to support inputs like "`A.B>>C`".
        if len(text) > 2 and text[0] == text[-1] and text[0] in {"`", '"', "'"}:
            text = text[1:-1].strip()

        match = REACTION_SMILES.search(text)
        if match:
            return match.group(1).strip()

        # Fallback: accept a whitespace-free token containing ">>".
        if " " not in text and "\t" not in text and "\n" not in text:
            return text
        return None

    def _is_conditions_intent(self, content: str, reaction_smiles: str) -> bool:
        normalized = (content or "").strip()
        if not normalized:
            return False

        # If user only pasted a reaction SMILES, assume they want conditions.
        wrapper_stripped = normalized
        if (
            len(wrapper_stripped) > 2
            and wrapper_stripped[0] == wrapper_stripped[-1]
            and wrapper_stripped[0] in {"`", '"', "'"}
        ):
            wrapper_stripped = wrapper_stripped[1:-1].strip()
        if wrapper_stripped == reaction_smiles:
            return True

        lowered = normalized.lower()
        return any(keyword in lowered for keyword in CONDITION_KEYWORDS)

    def _handle_conditions_command(self, raw_command: str, args: List[str]) -> None:
        if not args:
            self.append_chat(
                "System",
                "Usage: `/conditions <reaction_smiles> [optional constraint text]`.\n"
                "Tip: you can also paste a bare reaction SMILES to trigger this mode.",
            )
            return

        joined = " ".join(args).strip()
        reaction_smiles = self._extract_reaction_smiles(joined)
        if not reaction_smiles:
            self.append_chat(
                "System",
                "No reaction SMILES detected. Expected `reactants>>products`.\n"
                "Example: `/conditions Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1`",
            )
            return

        inline_text = joined.replace(reaction_smiles, " ").strip() or None
        inline_spec = build_constraint_spec(text=inline_text) if inline_text else ConstraintSpec()
        active_spec = merge_specs(self.constraint_spec, inline_spec)

        self.append_chat("You", raw_command)
        self._start_conditions_workflow(reaction_smiles, active_spec, inline_text=inline_text)

    def _maybe_start_conditions_workflow(self, content: str) -> bool:
        reaction_smiles = self._extract_reaction_smiles(content)
        if not reaction_smiles:
            return False
        if not self._is_conditions_intent(content, reaction_smiles):
            return False

        inline_text = content.replace(reaction_smiles, " ").strip() or None
        inline_spec = build_constraint_spec(text=inline_text) if inline_text else ConstraintSpec()
        active_spec = merge_specs(self.constraint_spec, inline_spec)

        self.append_chat("You", content)
        self._start_conditions_workflow(reaction_smiles, active_spec, inline_text=inline_text)
        return True

    def _start_conditions_workflow(
        self,
        reaction_smiles: str,
        constraints: ConstraintSpec,
        *,
        inline_text: Optional[str] = None,
    ) -> None:
        """Run a comprehensive, multi-database conditions workflow in the background."""
        constraint_summary = constraints.formatted_summary()
        constraint_display = constraint_summary if constraint_summary != "none" else "None"
        self.append_chat(
            "System",
            "Auto-conditions mode started.\n"
            f"- reaction: {reaction_smiles}\n"
            f"- constraints: {constraint_display}",
        )

        # Best-effort image preview for the query reaction.
        self.render_image_from_smiles(
            "reaction",
            reaction_smiles,
            source="query",
            silent=True,
        )

        self.set_input_enabled(False)
        self.start_spinner("Searching rules/ML/protocol/HTE")

        worker = Worker(
            self._run_conditions_workflow,
            reaction_smiles,
            constraints,
            inline_text=inline_text,
        )
        worker.signals.progress.connect(self.update_spinner_message)
        worker.signals.result.connect(self._handle_conditions_workflow_result)
        worker.signals.error.connect(self.handle_agent_error)
        worker.signals.finished.connect(self.reset_after_task)
        self.thread_pool.start(worker)

    def update_spinner_message(self, message: str) -> None:
        """Update the message shown next to the spinner."""
        self.spinner_message = message
        self._update_spinner()

    def _handle_conditions_workflow_result(self, payload: Dict[str, Any]) -> None:
        response = str(payload.get("response_text") or "").strip()
        details_path = payload.get("details_path")
        if details_path:
            self.last_details_path = details_path
            self.view_results_btn.setEnabled(True)
            response = response + f"\n\nDetails JSON saved to: {details_path}"

        if not response:
            self.append_chat("System", "Auto-conditions workflow returned no output.")
            return

        self.append_chat("Agent", response)
        self.append_log("Auto-conditions report ready.")

    def _run_conditions_workflow(
        self,
        reaction_smiles: str,
        constraints: ConstraintSpec,
        *,
        inline_text: Optional[str] = None,
        progress_callback: Optional[callable] = None,
    ) -> Dict[str, Any]:
        """Blocking worker: collect evidence from tools and synthesize a report."""
        from langchain_core.messages import HumanMessage, SystemMessage

        def report_progress(msg: str) -> None:
            if progress_callback:
                progress_callback(msg)

        from chem_assistant.chemtools_agent import get_llm_client
        from chem_assistant.chemtools_wrapper import (
            detect_reaction_family_tool,
            hte_recommend_tool,
            protocol_recommendation_tool,
            recommend_conditions_tool,
            rule_based_conditions_tool,
            unified_recommender_tool,
        )

        def safe_invoke(tool_obj: Any, tool_input: Dict[str, Any], label: str) -> Dict[str, Any]:
            report_progress(f"Running {label}...")
            try:
                result = tool_obj.invoke(tool_input)
            except Exception as exc:
                return {"success": False, "error": str(exc)}
            return result if isinstance(result, dict) else {"success": True, "data": result}

        constraint_payload: Dict[str, Any] = {}
        # ... (rest of the logic)
        if constraints.allow_metals:
            constraint_payload["allow_metals"] = sorted(constraints.allow_metals)
        if constraints.exclude_metals:
            constraint_payload["exclude_metals"] = sorted(constraints.exclude_metals)
        if constraints.prefer_metals:
            constraint_payload["prefer_metals"] = sorted(constraints.prefer_metals)
        if constraints.search_all_families:
            constraint_payload["search_all_families"] = True
        if constraints.constraint_rules:
            constraint_payload["constraint_rules"] = dict(constraints.constraint_rules)

        # Query all sources (best-effort; each call may fail independently).
        family = safe_invoke(detect_reaction_family_tool, {"reaction_smiles": reaction_smiles}, "reaction family detection")

        rule = safe_invoke(
            rule_based_conditions_tool,
            {"reaction_smiles": reaction_smiles, "include_summary": True},
            "rule-based engine"
        )

        ml = safe_invoke(
            recommend_conditions_tool,
            {
                "reaction_smiles": reaction_smiles,
                "k": 25,
                "max_variants": 3,
                "rerank_strategy": "rule",
                **constraint_payload,
            },
            "ML/DRFP recommender"
        )

        unified = safe_invoke(
            unified_recommender_tool,
            {
                "reaction_smiles": reaction_smiles,
                "top_k": 5,
                "min_similarity": 0.3,
                "validate_rules": True,
            },
            "unified search"
        )

        protocols = safe_invoke(
            protocol_recommendation_tool,
            {
                "reaction_smiles": reaction_smiles,
                "k": 3,
                "min_similarity": 0.3,
                "use_smarts_filter": True,
            },
            "protocol database"
        )

        # HTE: pull a larger set then filter/reorder by metal constraints.
        hte_raw = safe_invoke(
            hte_recommend_tool,
            {
                "reaction_smiles": reaction_smiles,
                "top_k": 20,
                "min_experiments": 1,
            },
            "HTE database"
        )
        hte = self._apply_hte_constraints(hte_raw, constraints, top_k=5)

        report_progress("Synthesizing report")

        details = {
            "reaction_smiles": reaction_smiles,
            "inline_text": inline_text,
            "constraints": {
                "allow_metals": sorted(constraints.allow_metals),
                "exclude_metals": sorted(constraints.exclude_metals),
                "prefer_metals": sorted(constraints.prefer_metals),
                "search_all_families": constraints.search_all_families,
                "constraint_rules": dict(constraints.constraint_rules),
            },
            "results": {
                "family": family,
                "rule": rule,
                "ml": ml,
                "unified": unified,
                "protocols": protocols,
                "hte": hte,
            },
        }

        details_jsonable = self._to_jsonable(details)
        details_path = (
            Path(tempfile.gettempdir())
            / f"chemtools_conditions_{uuid4().hex}.json"
        )
        details_path.write_text(json.dumps(details_jsonable, indent=2), encoding="utf-8")

        summary = self._build_conditions_report(
            reaction_smiles=reaction_smiles,
            constraints=constraints,
            family=family,
            rule=rule,
            ml=ml,
            protocols=protocols,
            unified=unified,
            hte=hte,
        )

        llm_report = self._maybe_llm_summarize_conditions(
            reaction_smiles=reaction_smiles,
            constraints=constraints,
            family=family,
            rule=rule,
            ml=ml,
            protocols=protocols,
            unified=unified,
            hte=hte,
            get_llm_client=get_llm_client,
            SystemMessage=SystemMessage,
            HumanMessage=HumanMessage,
        )

        return {
            "response_text": llm_report or summary,
            "details_path": str(details_path),
        }

    def _apply_hte_constraints(
        self,
        payload: Dict[str, Any],
        constraints: ConstraintSpec,
        *,
        top_k: int = 5,
    ) -> Dict[str, Any]:
        """Filter/reorder HTE recommendations by allow/exclude/prefer metals."""
        if not isinstance(payload, dict):
            return payload
        if not payload.get("success"):
            return payload
        recs = payload.get("recommendations")
        if not isinstance(recs, list) or not recs:
            return payload
        if not (constraints.allow_metals or constraints.exclude_metals or constraints.prefer_metals):
            payload["recommendations"] = self._to_jsonable(recs[:top_k])
            return payload

        def matches_any(catalyst: str, metals: set[str]) -> bool:
            upper = (catalyst or "").upper()
            return any(metal in upper for metal in metals)

        allow = set(constraints.allow_metals or set())
        exclude = set(constraints.exclude_metals or set())
        prefer = set(constraints.prefer_metals or set())

        filtered: List[Dict[str, Any]] = []
        for rec in recs:
            if not isinstance(rec, dict):
                continue
            catalyst = str(rec.get("catalyst") or "")
            if exclude and matches_any(catalyst, exclude):
                continue
            if allow and not matches_any(catalyst, allow):
                continue
            filtered.append(rec)

        if not filtered:
            filtered = [rec for rec in recs if isinstance(rec, dict)]

        if prefer:
            preferred = [rec for rec in filtered if matches_any(str(rec.get("catalyst") or ""), prefer)]
            others = [rec for rec in filtered if rec not in preferred]
            filtered = preferred + others

        payload["recommendations"] = self._to_jsonable(filtered[:top_k])
        payload["hte_constraints_applied"] = constraints.formatted_summary()
        return payload

    def _build_conditions_report(
        self,
        *,
        reaction_smiles: str,
        constraints: ConstraintSpec,
        family: Dict[str, Any],
        rule: Dict[str, Any],
        ml: Dict[str, Any],
        protocols: Dict[str, Any],
        unified: Dict[str, Any],
        hte: Dict[str, Any],
    ) -> str:
        """Create a concise, human-readable summary from collected tool outputs."""
        lines: List[str] = []
        lines.append("COMPREHENSIVE CONDITIONS REPORT")
        lines.append("=" * 72)
        lines.append(f"Reaction: {reaction_smiles}")
        summary = constraints.formatted_summary()
        if summary != "none":
            lines.append(f"Constraints: {summary}")
        lines.append("")

        # Family
        if family.get("success"):
            fam = family.get("family")
            conf = family.get("confidence")
            method = family.get("method")
            lines.append(f"Detected family: {fam} (confidence={conf}, method={method})")
        else:
            lines.append(f"Detected family: unavailable ({family.get('error')})")
        lines.append("")

        # Rule
        lines.append("[Rule DB]")
        if rule.get("success"):
            rule_summary = str(rule.get("summary") or "").strip()
            if rule_summary:
                lines.append(rule_summary)
            else:
                lines.append("Rule engine returned no summary.")
        else:
            lines.append(f"Rule engine error: {rule.get('error')}")
        lines.append("")

        # ML / DRFP
        lines.append("[ML/DRFP Precedents]")
        if ml.get("success"):
            rec = ml.get("recommendation") or {}
            core = rec.get("core")
            base = rec.get("base") or rec.get("base_uid")
            solv = rec.get("solvent") or rec.get("solvent_uid")
            temp = rec.get("T_C")
            time_h = rec.get("time_h")
            conf = rec.get("confidence")
            lines.append(f"- core: {core}")
            lines.append(f"- base: {base}")
            lines.append(f"- solvent: {solv}")
            if temp is not None:
                lines.append(f"- temperature_C: {temp}")
            if time_h is not None:
                lines.append(f"- time_h: {time_h}")
            if conf is not None:
                lines.append(f"- confidence: {conf}")
            alt = ml.get("alternatives") or {}
            cores = alt.get("cores") or []
            if cores:
                lines.append(f"- alternative cores: {', '.join(str(c[0]) for c in cores[:5] if c)}")
        else:
            lines.append(f"ML/DRFP recommendation error: {ml.get('error')}")
        lines.append("")

        # HTE
        lines.append("[HTE (66k experiments)]")
        if hte.get("success") and (hte.get("matching_experiments") or 0) > 0:
            lines.append(
                f"Matches: {hte.get('matching_experiments')} | predicted type: {hte.get('predicted_reaction_type')}"
            )
            recs = hte.get("recommendations") or []
            if isinstance(recs, list) and recs:
                lines.append("Top HTE conditions:")
                for rec in recs[:5]:
                    if not isinstance(rec, dict):
                        continue
                    lines.append(
                        f"- #{rec.get('rank')}: {rec.get('catalyst')} / {rec.get('ligand')} | "
                        f"{rec.get('base')} | {rec.get('solvent')} | "
                        f"succ={rec.get('success_rate')}% (n={rec.get('num_experiments')}) | "
                        f"avg_yield={rec.get('avg_yield')} | z={rec.get('avg_z_score')}"
                    )
            else:
                lines.append("No HTE condition list returned.")
        elif hte.get("success"):
            lines.append("No matching HTE experiments found for this substrate pairing.")
        else:
            lines.append(f"HTE error: {hte.get('error')}")
        lines.append("")

        # Protocol DB
        lines.append("[Protocol DB]")
        if protocols.get("success"):
            count = protocols.get("count") or len(protocols.get("recommendations") or [])
            lines.append(f"Protocols found: {count}")
        else:
            hint = protocols.get("hint")
            lines.append(f"Protocol search unavailable: {protocols.get('error')}")
            if hint:
                lines.append(f"Hint: {hint}")
        lines.append("")

        # Unified (rules/protocols search)
        lines.append("[Unified DRFP Search]")
        if unified.get("success"):
            recs = unified.get("recommendations") or []
            if isinstance(recs, list) and recs:
                for item in recs[:5]:
                    if not isinstance(item, dict):
                        continue
                    lines.append(
                        f"- #{item.get('rank')}: {item.get('name')} ({item.get('source_type')}) "
                        f"sim={item.get('similarity')}"
                    )
            else:
                lines.append("Unified search returned no items.")
        else:
            lines.append(f"Unified search error: {unified.get('error')}")
        lines.append("")

        # Simple deterministic pick
        lines.append("[Recommended Starting Point]")
        best = None
        if hte.get("success") and isinstance(hte.get("recommendations"), list) and hte["recommendations"]:
            best = hte["recommendations"][0]
            lines.append(
                f"HTE-backed: {best.get('catalyst')} / {best.get('ligand')} | "
                f"{best.get('base')} | {best.get('solvent')} "
                f"(succ={best.get('success_rate')}%, n={best.get('num_experiments')})"
            )
        elif rule.get("success") and rule.get("summary"):
            lines.append("Rule-backed: use the Rule DB recommendation above as the primary starting point.")
        elif ml.get("success") and ml.get("recommendation"):
            rec = ml.get("recommendation") or {}
            lines.append(
                f"Precedent-backed: core={rec.get('core')}, base={rec.get('base') or rec.get('base_uid')}, "
                f"solvent={rec.get('solvent') or rec.get('solvent_uid')} (confidence={rec.get('confidence')})"
            )
        else:
            lines.append("No successful recommendation sources available; check tool errors above.")

        return "\n".join(str(line) for line in lines if line is not None)

    def _maybe_llm_summarize_conditions(
        self,
        *,
        reaction_smiles: str,
        constraints: ConstraintSpec,
        family: Dict[str, Any],
        rule: Dict[str, Any],
        ml: Dict[str, Any],
        protocols: Dict[str, Any],
        unified: Dict[str, Any],
        hte: Dict[str, Any],
        get_llm_client: Any,
        SystemMessage: Any,
        HumanMessage: Any,
    ) -> Optional[str]:
        """Use the configured LLM (no tools) to synthesize a ranked recommendation."""
        try:
            llm = get_llm_client(self.llm_provider, self.llm_model, self.llm_temperature)
        except Exception:
            return None

        # Build a compact evidence payload (avoid giant vector dumps).
        evidence: Dict[str, Any] = {
            "reaction_smiles": reaction_smiles,
            "constraints": constraints.formatted_summary(),
            "family": family if family.get("success") else {"success": False, "error": family.get("error")},
            "rule_summary": rule.get("summary") if rule.get("success") else {"error": rule.get("error")},
            "ml_recommendation": ml.get("recommendation") if ml.get("success") else {"error": ml.get("error")},
            "ml_alternatives": (ml.get("alternatives") or {}) if ml.get("success") else None,
            "hte": {
                "predicted_reaction_type": hte.get("predicted_reaction_type"),
                "matching_experiments": hte.get("matching_experiments"),
                "recommendations": hte.get("recommendations"),
                "error": None if hte.get("success") else hte.get("error"),
            },
            "protocols_status": {
                "success": bool(protocols.get("success")),
                "error": protocols.get("error"),
                "hint": protocols.get("hint"),
            },
            "unified_top": unified.get("recommendations") if unified.get("success") else {"error": unified.get("error")},
        }

        system = (
            "You are an expert synthetic chemist. You will be given tool outputs from:\n"
            "- Rule DB (deterministic rules)\n"
            "- ML/DRFP precedent recommender\n"
            "- Protocol DB (may be unavailable)\n"
            "- HTE database (66k experiments)\n"
            "- A deterministic planner fusion snapshot\n\n"
            "Task:\n"
            "1) Briefly summarize the reaction family and constraints.\n"
            "2) Compare evidence across sources.\n"
            "3) Provide 2-3 ranked, executable condition recommendations.\n"
            "4) For each recommendation: catalyst/ligand/base/solvent, temperature/time if available, and cite which sources support it.\n"
            "5) Do not invent missing values; if not available, say 'not provided'.\n"
        )
        user = "Tool evidence (JSON):\n" + json.dumps(self._to_jsonable(evidence), indent=2)

        try:
            response = llm.invoke([SystemMessage(content=system), HumanMessage(content=user)])
        except Exception:
            return None

        content = getattr(response, "content", None)
        return str(content).strip() if content else None

    def _to_jsonable(self, obj: Any) -> Any:
        """Best-effort conversion of tool outputs into JSON-serializable structures."""
        if obj is None:
            return None
        if isinstance(obj, (str, int, float, bool)):
            return obj
        if isinstance(obj, Path):
            return str(obj)
        if isinstance(obj, dict):
            return {str(k): self._to_jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple, set)):
            return [self._to_jsonable(v) for v in obj]
        if hasattr(obj, "model_dump"):
            try:
                return self._to_jsonable(obj.model_dump())  # type: ignore[attr-defined]
            except Exception:
                return str(obj)
        if hasattr(obj, "item"):
            try:
                return obj.item()  # type: ignore[no-any-return]
            except Exception:
                pass
        try:
            json.dumps(obj)
            return obj
        except TypeError:
            return str(obj)

    def _ensure_agent(self) -> bool:
        if self.agent is not None:
            return True
        try:
            self.agent = ChemToolsAgent(
                verbose=True,
                provider=self.llm_provider,
                model=self.llm_model,
                temperature=self.llm_temperature,
                proactive=self.proactive_enabled,
                proactive_top_k=self.proactive_top_k,
                proactive_max_protocols=self.proactive_max_protocols,
                proactive_build_protocols=self.proactive_build_protocols,
            )
        except Exception as exc:  # pragma: no cover - UI feedback
            self.append_chat(
                "System",
                "LLM agent is not available.\n"
                f"Reason: {exc}\n\n"
                "You can still use local commands:\n"
                "- `/taxonomy` (show taxonomy status)\n"
                "- `/terms <SMILES>` (evaluate chemistry terms)\n"
                "- `/conditions <reaction_smiles>` (multi-tool conditions search)\n"
                "- `image <reaction|molecule> <SMILES>`",
            )
            return False
        return True

    def view_last_results(self) -> None:
        if not self.last_details_path or not Path(self.last_details_path).exists():
            QMessageBox.warning(self, "No results", "No recent condition results found.")
            return
        
        try:
            data = json.loads(Path(self.last_details_path).read_text(encoding="utf-8"))
            dialog = JsonViewerDialog(data, title="Condition Search Details", parent=self)
            dialog.exec()
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Failed to load results: {exc}")

    def configure_llm_mode(self) -> None:
        dialog = LLMConfigDialog(
            provider=self.llm_provider,
            model=self.llm_model,
            temperature=self.llm_temperature,
            proactive=self.proactive_enabled,
            proactive_top_k=self.proactive_top_k,
            proactive_max_protocols=self.proactive_max_protocols,
            proactive_build_protocols=self.proactive_build_protocols,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        config = dialog.get_config()
        self.llm_provider = config.provider
        self.llm_model = config.model
        self.llm_temperature = config.temperature
        self.proactive_enabled = config.proactive
        self.proactive_top_k = config.proactive_top_k
        self.proactive_max_protocols = config.proactive_max_protocols
        self.proactive_build_protocols = config.proactive_build_protocols

        # Recreate agent lazily with the new settings on next send.
        self.agent = None
        proactive_status = "on" if self.proactive_enabled else "off"
        self.append_log(
            "LLM mode set: "
            f"{self.llm_provider}/{self.llm_model} (T={self.llm_temperature:g}); "
            f"proactive preflight {proactive_status}"
        )

    def show_taxonomy_status(self) -> None:
        self._handle_taxonomy_command()

    def _handle_taxonomy_command(self) -> None:
        try:
            from chemtools.taxonomy.reaction_catalog import load_reaction_catalog
            from chemtools.taxonomy.reagent_v2 import ReagentTaxonomyV2

            rxn_defs, _ = load_reaction_catalog()
            reagent_tax = ReagentTaxonomyV2.from_path()

            rxn_count = len(rxn_defs)
            reagent_count = len(list(reagent_tax.iter_families()))

            self.append_chat(
                "System",
                "Taxonomy v2 loaded.\n"
                f"- reaction_types: {rxn_count}\n"
                f"- reagent_families: {reagent_count}",
            )
        except Exception as exc:  # pragma: no cover - UI feedback
            self.append_chat("System", f"Taxonomy v2 load failed: {exc}")

    def _handle_terms_command(self, cmd: str, args: List[str]) -> None:
        if cmd == "term":
            if len(args) < 2:
                self.append_chat("System", "Usage: `/term <term_id> <SMILES>`")
                return
            term_ids = [args[0]]
            smiles = " ".join(args[1:]).strip()
        else:
            if not args:
                self.append_chat("System", "Usage: `/terms <SMILES>`")
                return
            term_ids = None
            smiles = " ".join(args).strip()

        try:
            from chemtools.util.functional_groups import detect_all, summarize_functional_groups

            if term_ids:
                results = detect_all(smiles)
                hits = [tid for tid in term_ids if results.get(tid)]
                if not hits:
                    self.append_chat("System", f"No matches for requested terms: {', '.join(term_ids)}")
                    return
                lines = [f"- {tid}" for tid in sorted(hits)]
                self.append_chat("System", "Matches:\n" + "\n".join(lines))
            else:
                summary = summarize_functional_groups(smiles)
                self.append_chat("System", f"Functional Groups:\n{summary}")
        except Exception as exc:  # pragma: no cover - UI feedback
            self.append_chat("System", f"Term evaluation failed: {exc}")

    def handle_image_command(self, command: str) -> None:
        parsed = self._extract_image_request(command)
        if not parsed:
            self.append_chat(
                "System",
                "Usage: image <reaction|molecule> <SMILES>. Example: "
                "image reaction BrC1(c2ccccc2)CC1.c1ccc(B(O)O)cc1>>...",
            )
            return
        target, smiles = parsed
        self._execute_image_request(target, smiles, log_to_chat=True)

    def _execute_image_request(
        self,
        target: str,
        smiles: str,
        *,
        log_to_chat: bool = False,
    ) -> None:
        if not smiles:
            self.append_chat("System", "Provide a SMILES string to render.")
            return
        if self.render_image_from_smiles(target, smiles, source="command"):
            if log_to_chat:
                self.append_chat(
                    "System",
                    f"{target.title()} image opened in a pop-out preview.",
                )

    def handle_agent_result(self, payload: Tuple[str, List[BaseMessage]]) -> None:
        try:
            response, history = payload
        except ValueError:
            self.append_log("Unexpected agent payload.")
            return
        self.history = history
        clean_response = self._process_agent_image_markup(response)
        
        # Auto-detect and render SMILES from response
        self._auto_render_from_response(response)
        
        self.append_chat("Agent", clean_response)
        self.append_log("Agent response received.")

    def _auto_render_from_response(self, text: str) -> None:
        """Automatically detect and render SMILES strings from agent response."""
        if not text:
            return
        
        # Pattern to match SMILES strings in various contexts
        # Look for reaction SMILES (contains >>)
        reaction_patterns = [
            r'(?:reaction|rxn)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[\S]+)',
            r'(?:SMILES|smiles)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[\S]+)',
            r'`([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)`',
            r'"([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)"',
            r'\b([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{10,}>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)\b',
        ]
        
        # Look for molecule SMILES (no >>)
        molecule_patterns = [
            r'(?:compound|molecule|structure)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.]{3,})',
            r'(?:SMILES|smiles)[\s:=]+([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{2,})',
            r'`([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{5,})`',
        ]
        
        # Try to find reaction SMILES first (higher priority)
        for pattern in reaction_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                smiles = match.group(1).strip()
                # Validate it looks like a reaction
                if '>>' in smiles and len(smiles) > 10:
                    self.render_image_from_smiles(
                        "reaction",
                        smiles,
                        source="auto-detected",
                        silent=True
                    )
                    return
        
        # If no reaction found, try molecule SMILES
        for pattern in molecule_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                smiles = match.group(1).strip()
                # Skip if it contains reaction arrow or looks like plain text
                if '>>' not in smiles and len(smiles) >= 5:
                    # Basic validation: should start with letter or bracket
                    # and contain typical SMILES characters
                    if (smiles[0].isalpha() or smiles[0] in '[(') and \
                       any(c in smiles for c in '()[]=#123456'):
                        self.render_image_from_smiles(
                            "molecule",
                            smiles,
                            source="auto-detected",
                            silent=True
                        )
                        return

    def _process_agent_image_markup(self, text: str) -> str:
        if not text:
            return text

        def directive_replacer(match: re.Match[str]) -> str:
            kind = match.group(1).lower()
            smiles = match.group(2).strip()
            normalized_kind = "reaction" if "reaction" in kind else "molecule"
            success = self.render_image_from_smiles(
                normalized_kind,
                smiles,
                source="agent",
                silent=True,
            )
            return "[image rendered]" if success else "[image unavailable]"

        interim = IMAGE_MARKUP.sub(directive_replacer, text)

        def markdown_replacer(match: re.Match[str]) -> str:
            url = match.group(1)
            if self._handle_markdown_image_url(url):
                return "[image rendered]"
            return match.group(0)

        interim = MARKDOWN_IMAGE.sub(markdown_replacer, interim)
        return interim

    def _extract_image_request(self, content: str) -> Optional[Tuple[str, str]]:
        normalized = content.strip()
        if normalized.startswith("/"):
            normalized = normalized[1:].lstrip()
        for pattern in IMAGE_COMMAND_PATTERNS:
            match = re.match(pattern, normalized, re.IGNORECASE)
            if not match:
                continue
            type_token = match.group(1) or ""
            smiles = match.group(2).strip() if match.lastindex and match.lastindex >= 2 else ""
            target = self._normalize_image_target(type_token, smiles)
            if not smiles:
                return None
            if target is None:
                self.append_chat(
                    "System",
                    "Unknown image target. Use 'reaction' or 'molecule'.",
                )
                return None
            return target, smiles
        return None

    def _normalize_image_target(
        self,
        type_token: str,
        smiles: str,
    ) -> Optional[str]:
        token = (type_token or "").strip().lower()
        if token in {"reaction", "rxn"}:
            return "reaction"
        if token in {"molecule", "compound", "mol"}:
            return "molecule"
        if not token:
            return "reaction" if ">" in smiles else "molecule"
        return None

    def _handle_markdown_image_url(self, url: str) -> bool:
        try:
            parsed = urlparse(url)
        except Exception:
            return False
        query = parse_qs(parsed.query)
        smiles = ""
        for key in URL_SMILES_KEYS:
            values = query.get(key)
            if values:
                smiles = values[0].strip()
                if smiles:
                    break
        if not smiles:
            return False
        if not smiles:
            return False
        target = "reaction" if ">" in smiles else "molecule"
        return self.render_image_from_smiles(
            target,
            smiles,
            source="agent",
            silent=True,
        )

    def handle_agent_error(self, trace: str) -> None:
        QMessageBox.critical(self, "Agent error", trace)
        self.append_log(trace)

    def reset_after_task(self) -> None:
        self.set_input_enabled(True)
        self.statusBar().showMessage("Ready")
        self.stop_spinner()

    def set_input_enabled(self, enabled: bool) -> None:
        self.send_btn.setEnabled(enabled)
        self.input_edit.setEnabled(enabled)

    def clear_chat(self) -> None:
        self.chat_view.clear()
        self.history = []
        self.append_log("Chat cleared.")

    def render_image_from_smiles(
        self,
        target: str,
        smiles: str,
        *,
        source: str,
        silent: bool = False,
    ) -> bool:
        try:
            image_path = self._generate_image_file(target, smiles)
        except ValueError as exc:
            if not silent:
                QMessageBox.warning(self, "Invalid SMILES", str(exc))
            else:
                self.append_log(f"Image render failed: {exc}")
            return False
        except RuntimeError as exc:
            if not silent:
                QMessageBox.warning(self, "Rendering unavailable", str(exc))
            else:
                self.append_log(str(exc))
            return False
        except Exception as exc:
            if not silent:
                QMessageBox.warning(self, "Image rendering error", str(exc))
            else:
                self.append_log(f"Image render error: {exc}")
            return False

        try:
            self.display_image(image_path, f"{target.title()} ({source})")
        except ValueError as exc:
            if not silent:
                QMessageBox.warning(self, "Image preview error", str(exc))
            else:
                self.append_log(f"Image preview error: {exc}")
            return False
        return True

    def _generate_image_file(self, target: str, smiles: str) -> Path:
        if not smiles:
            raise ValueError("SMILES string cannot be empty.")
        destination = Path(tempfile.gettempdir()) / f"chemtools_{target}_{uuid4().hex}.png"
        if target == "reaction":
            render_reaction_image(smiles, destination, image_format="png")
        else:
            render_molecule_image(smiles, destination, image_format="png")
        return destination

    def display_image(self, image_path: Path, caption: str) -> None:
        pixmap = QPixmap(str(image_path))
        if pixmap.isNull():
            raise ValueError("Unable to load rendered image.")
        preview = ImagePreviewWindow(pixmap, caption, parent=self)
        self._register_image_window(preview)
        preview.show()
        self.current_image_path = image_path
        self.append_log(f"Image preview opened: {caption or image_path.name}")

    def _register_image_window(self, window: ImagePreviewWindow) -> None:
        self._image_windows.append(window)

        def _cleanup(_: QObject, target: ImagePreviewWindow = window) -> None:
            self._unregister_image_window(target)

        window.destroyed.connect(_cleanup)

    def _unregister_image_window(self, window: ImagePreviewWindow) -> None:
        try:
            self._image_windows.remove(window)
        except ValueError:
            pass

    # ------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------ #

    def manage_constraints(self) -> None:
        dialog = ConstraintDialog(self.constraint_spec, self.constraint_text, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:  # type: ignore[attr-defined]
            self.constraint_spec = dialog.build_spec()
            self.constraint_text = dialog.get_constraint_text()
            self.append_log("Constraints updated.")
            self._update_constraint_summary()

    def _update_constraint_summary(self) -> None:
        summary = self.constraint_spec.formatted_summary()
        display = summary or ""
        if not display or display.strip().lower() == "none":
            display = "None"
        self.append_log(f"Constraints: {display}")

        prompt = format_constraints_for_prompt(self.constraint_spec)
        if prompt:
            self.append_log(f"Constraint prompt: {prompt}")

    # ------------------------------------------------------------------ #
    # Spinner helpers
    # ------------------------------------------------------------------ #

    def start_spinner(self, message: str) -> None:
        if not hasattr(self, "spinner_chars"):
            self.spinner_chars = ["|", "/", "-", "\\"]
            self.spinner_index = 0
            self.spinner_timer = QTimer(self)
            self.spinner_timer.timeout.connect(self._update_spinner)
        self.spinner_message = message
        self.statusBar().showMessage(f"{message}... {self.spinner_chars[self.spinner_index]}")
        self.spinner_timer.start(120)

    def stop_spinner(self) -> None:
        if hasattr(self, "spinner_timer"):
            self.spinner_timer.stop()
        self.statusBar().showMessage("Ready")

    def _update_spinner(self) -> None:
        self.spinner_index = (self.spinner_index + 1) % len(self.spinner_chars)
        self.statusBar().showMessage(
            f"{self.spinner_message}... {self.spinner_chars[self.spinner_index]}"
        )

    # ------------------------------------------------------------------ #
    # Drag and Drop
    # ------------------------------------------------------------------ #

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasFormat("text/plain") or event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:
        text = ""
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                path = Path(url.toLocalFile())
                if path.suffix.lower() in {".smiles", ".txt", ".mol", ".sdf"}:
                    try:
                        text = path.read_text(encoding="utf-8").strip()
                        break
                    except Exception:
                        pass
        elif event.mimeData().hasFormat("text/plain"):
            text = event.mimeData().text().strip()

        if text:
            self.input_edit.setPlainText(text)
            self.append_log(f"Dropped content loaded into input.")
            event.acceptProposedAction()
