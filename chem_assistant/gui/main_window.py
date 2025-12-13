"""Main Qt window for the ChemTools assistant."""

from __future__ import annotations

import html
import json
import re
import shlex
import tempfile
import traceback
from pathlib import Path
from typing import List, Optional, Tuple
from urllib.parse import parse_qs, urlparse
from uuid import uuid4

from PyQt6.QtCore import QObject, QRunnable, Qt, QThreadPool, pyqtSignal, QTimer
from PyQt6.QtWidgets import (
    QDialog,
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
from chem_assistant.chemtools_wrapper import (
    get_tool_descriptions,
    clear_recommendation_cache,
    recommendation_cache_stats,
)
from chem_assistant.constraint_parser import (
    ConstraintSpec,
    format_constraints_for_prompt,
)
from chem_assistant.gui.dialogs import (
    ConstraintDialog,
    LLMConfigDialog,
    RuleBuilderAutofillDialog,
    RuleBuilderDialog,
    ProtocolDraftDialog,
)
from chemtools.visualization import render_molecule_image, render_reaction_image

IMAGE_MARKUP = re.compile(
    r"\[\[(reaction|molecule)_image:(.+?)\]\]", re.IGNORECASE | re.DOTALL
)
MARKDOWN_IMAGE = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")
PROTOCOL_PROMPT = re.compile(
    r"\[\[protocol_draft_request(?::(.*?))?\]\]", re.IGNORECASE | re.DOTALL
)
IMAGE_COMMAND_PATTERNS = [
    r"^(?:/)?image\s+(reaction|molecule|compound)\s*[:=]\s*(.+)$",
    r"^(?:/)?image\s+(?:for\s+)?(reaction|molecule|compound)\s*[:=]\s*(.+)$",
    r"^(?:show|display)\s+(?:the\s+)?image(?:\s+of|\s+for)?\s*(?:a|the)?\s*(reaction|molecule|compound)?[:=]?\s*(.+)$",
    r"^(?:show|display)\s+(?:me\s+)?(?:an?\s+|the\s+)?image\s*(?:of|for)?\s*(reaction|molecule|compound)?[:=]?\s*(.+)$",
]
URL_SMILES_KEYS = ("smiles", "model", "structure", "mol", "target")


class WorkerSignals(QObject):
    """Signals for background tasks."""

    finished = pyqtSignal()
    error = pyqtSignal(str)
    result = pyqtSignal(object)


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
        self.history: List[BaseMessage] = []
        self.constraint_spec = ConstraintSpec()
        self.constraint_text = ""
        self.current_image_path: Optional[Path] = None
        self._image_windows: List[ImagePreviewWindow] = []
        self._startup_message = (startup_message or "").strip()

        self.thread_pool = QThreadPool()
        self._build_ui()
        self._update_constraint_summary()
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
        self.input_edit.setPlaceholderText("Type a question or command...")
        self.input_edit.setMinimumHeight(60)
        self.input_edit.setMaximumHeight(140)
        self.input_edit.setStyleSheet(
            "QTextEdit { background-color: #30343b; color: #f5f6fa; border: 1px solid #555; }"
        )
        input_font = QFont(self.font())
        input_font.setPointSize(input_font.pointSize() + 1)
        self.input_edit.setFont(input_font)
        self.input_edit.sendRequested.connect(self.send_message)
        central_layout.addWidget(self.input_edit)

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

        cache_btn = QPushButton("Cache Status")
        cache_btn.clicked.connect(self.show_cache_stats)
        button_row.addWidget(cache_btn)

        llm_btn = QPushButton("LLM Mode")
        llm_btn.clicked.connect(self.configure_llm_mode)
        button_row.addWidget(llm_btn)

        cache_clear_btn = QPushButton("Clear Cache")
        cache_clear_btn.clicked.connect(self.clear_cache)
        button_row.addWidget(cache_clear_btn)

        tools_btn = QPushButton("List Tools")
        tools_btn.clicked.connect(self.show_tool_summary)
        button_row.addWidget(tools_btn)

        taxonomy_btn = QPushButton("Taxonomy Status")
        taxonomy_btn.clicked.connect(self.show_taxonomy_status)
        button_row.addWidget(taxonomy_btn)

        builder_btn = QPushButton("Rule Builder Editor")
        builder_btn.clicked.connect(self.open_rule_builder_dialog)
        button_row.addWidget(builder_btn)

        autofill_btn = QPushButton("Autofill Draft")
        autofill_btn.clicked.connect(self.open_autofill_dialog)
        button_row.addWidget(autofill_btn)

        protocol_btn = QPushButton("Protocol Draft")
        protocol_btn.clicked.connect(self.open_protocol_draft_dialog)
        button_row.addWidget(protocol_btn)

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

    def send_message(self) -> None:
        content = self.input_edit.toPlainText().strip()
        if not content:
            return
        self.input_edit.clear()
        if self.handle_local_command(content):
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

        parsed = self._extract_image_request(command_text)
        if parsed:
            target, smiles = parsed
            self._execute_image_request(target, smiles, log_to_chat=True)
            return True
        return False

    def _ensure_agent(self) -> bool:
        if self.agent is not None:
            return True
        try:
            self.agent = ChemToolsAgent(
                verbose=True,
                provider=self.llm_provider,
                model=self.llm_model,
                temperature=self.llm_temperature,
            )
        except Exception as exc:  # pragma: no cover - UI feedback
            self.append_chat(
                "System",
                "LLM agent is not available.\n"
                f"Reason: {exc}\n\n"
                "You can still use local commands:\n"
                "- `/taxonomy` (show taxonomy status)\n"
                "- `/terms <SMILES>` (evaluate chemistry terms)\n"
                "- `image <reaction|molecule> <SMILES>`",
            )
            return False
        return True

    def configure_llm_mode(self) -> None:
        dialog = LLMConfigDialog(
            provider=self.llm_provider,
            model=self.llm_model,
            temperature=self.llm_temperature,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        config = dialog.get_config()
        self.llm_provider = config.provider
        self.llm_model = config.model
        self.llm_temperature = config.temperature

        # Recreate agent lazily with the new settings on next send.
        self.agent = None
        self.append_log(
            f"LLM mode set: {self.llm_provider}/{self.llm_model} (T={self.llm_temperature:g})"
        )

    def show_taxonomy_status(self) -> None:
        self._handle_taxonomy_command()

    def _handle_taxonomy_command(self) -> None:
        try:
            from chemtools.taxonomy import load_registry

            registry = load_registry()
            term_count = len(list(registry.iter_chem_terms()))
            self.append_chat(
                "System",
                "Taxonomy loaded.\n"
                f"- taxonomy_version: {registry.manifest.taxonomy_version}\n"
                f"- schema_version: {registry.manifest.schema_version}\n"
                f"- chem_terms: {term_count}",
            )
        except Exception as exc:  # pragma: no cover - UI feedback
            self.append_chat("System", f"Taxonomy load failed: {exc}")

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
            from chemtools.taxonomy import load_registry
            from chemtools.taxonomy.terms import evaluate_terms_from_smiles

            registry = load_registry()
            results = evaluate_terms_from_smiles(smiles, term_ids=term_ids)
            hits = [term_id for term_id, ok in results.items() if ok]
            if not hits:
                self.append_chat("System", "No chem-term matches.")
                return

            lines = []
            for term_id in hits:
                term = registry.get_chem_term(term_id)
                if term is None:
                    lines.append(f"- {term_id}")
                else:
                    lines.append(f"- {term_id}: {term.name}")
            self.append_chat("System", "Chem-term matches:\n" + "\n".join(lines))
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
        return PROTOCOL_PROMPT.sub(self._handle_protocol_prompt_directive, interim)

    def _handle_protocol_prompt_directive(self, match: re.Match[str]) -> str:
        raw_payload = (match.group(1) or "").strip()
        data: Optional[dict] = None
        if raw_payload:
            try:
                data = json.loads(raw_payload)
            except json.JSONDecodeError:
                data = {"procedure_text": raw_payload}
        QTimer.singleShot(0, lambda: self.open_protocol_draft_dialog(data))
        return "[protocol draft dialog opened]"

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
    # Constraints / cache / tools
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

    def show_cache_stats(self) -> None:
        stats = recommendation_cache_stats()
        self.append_log(
            f"Cache status -> entries: {stats['entries']} | hits: {stats['hits']}"
        )
        QMessageBox.information(
            self,
            "Cache status",
            "\n".join(f"{k}: {v}" for k, v in stats.items()),
        )

    def clear_cache(self) -> None:
        clear_recommendation_cache()
        self.append_log("Recommendation cache cleared.")
        QMessageBox.information(self, "Cache", "Recommendation cache cleared.")

    def show_tool_summary(self) -> None:
        tools = get_tool_descriptions()
        lines = [f"{entry['name']}: {entry['description']}" for entry in tools]
        QMessageBox.information(self, "Available tools", "\n\n".join(lines))

    # ------------------------------------------------------------------ #
    # Builder dialogs
    # ------------------------------------------------------------------ #

    def open_rule_builder_dialog(self, initial_data: Optional[dict] = None) -> None:
        dialog = RuleBuilderDialog(self, initial_data=initial_data)
        dialog.exec()

    def open_autofill_dialog(self) -> None:
        dialog = RuleBuilderAutofillDialog(self)
        if dialog.exec() == QDialog.DialogCode.Accepted:  # type: ignore[attr-defined]
            if dialog.accepted_data:
                QMessageBox.information(
                    self,
                    "Autofill draft ready",
                    "Opening draft in editor...",
                )
                self.open_rule_builder_dialog(dialog.accepted_data)

    def open_protocol_draft_dialog(self, initial_data: Optional[dict] = None) -> None:
        dialog = ProtocolDraftDialog(self, initial_data=initial_data)
        if dialog.exec() != QDialog.DialogCode.Accepted:  # type: ignore[attr-defined]
            return
        payload = dialog.accepted_payload or {}
        draft = payload.get("draft") or {}
        metadata = draft.get("metadata") or {}
        name = metadata.get("name") or metadata.get("id") or "protocol draft"
        issues = payload.get("issues") or []
        summary = f"Protocol draft ready: {name}."
        if issues:
            summary += f" Outstanding issues: {len(issues)} (see preview for details)."
        self.append_chat("System", summary)
        if payload.get("llm_used"):
            meta = payload.get("llm_metadata") or {}
            self.append_chat(
                "System",
                "LLM extraction used "
                f"({meta.get('provider') or 'provider?'} / {meta.get('model') or 'model?'}).",
            )
        saved_path = payload.get("saved_path")
        if saved_path:
            self.append_chat("System", f"Draft saved to {saved_path}. Rebuild the protocol index to use it.")
        addition_sequence = payload.get("addition_sequence") or []
        if addition_sequence:
            first_steps = ", ".join(
                step.get("material", step.get("action", "step")).strip()
                for step in addition_sequence[:3]
                if step.get("material") or step.get("action")
            )
            if first_steps:
                self.append_chat("System", f"First addition steps: {first_steps}")

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
