"""Main Qt window for the ChemTools assistant."""

from __future__ import annotations

import html
import re
import shlex
import tempfile
import traceback
from pathlib import Path
from typing import Any, List, Optional, Tuple
from urllib.parse import parse_qs, urlparse
from uuid import uuid4

from PyQt6.QtCore import QObject, QRunnable, Qt, QThreadPool, pyqtSignal, QTimer, QSize
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
from chem_assistant.gui.dialogs import LLMConfigDialog
from chemtools.visualization import render_molecule_image, render_reaction_image
from chemtools.util.rdkit_helpers import parse_smiles
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

    _SCALE_FACTOR = 0.5

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
        self._max_width = max(1, int(pixmap.width() * self._SCALE_FACTOR))
        self._max_height = max(1, int(pixmap.height() * self._SCALE_FACTOR))

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
            max(420, self._max_width + 64),
            max(360, self._max_height + 96),
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
        bounded = QSize(
            min(target_size.width(), self._max_width),
            min(target_size.height(), self._max_height),
        )
        if not bounded.isValid():
            return
        scaled = self._original_pixmap.scaled(
            bounded,
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
        self.current_image_path: Optional[Path] = None
        self._image_windows: List[ImagePreviewWindow] = []
        self._startup_message = (startup_message or "").strip()

        self.thread_pool = QThreadPool()
        self._build_ui()
        self.setAcceptDrops(True)
        self.statusBar().showMessage(self._startup_message or "Ready")
        if self._startup_message:
            self.append_log(self._startup_message)

    def _apply_default_font(self) -> None:
        font = QFont(self.font())
        if font.pointSize() > 0:
            font.setPointSize(font.pointSize() + 2)
        else:
            font.setPointSize(12)
        self.setFont(font)

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

        self.input_edit = ChatInput()
        self.input_edit.setPlaceholderText(
            "Type a question, or paste SMILES. Use /taxonomy, /terms, or /image."
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

        llm_btn = QPushButton("LLM Mode")
        llm_btn.clicked.connect(self.configure_llm_mode)
        button_row.addWidget(llm_btn)

        button_row.addStretch()
        central_layout.addLayout(button_row)

        self.setCentralWidget(central)
        self.setStatusBar(QStatusBar())

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

        reaction_smiles = self._extract_reaction_smiles(content)
        if reaction_smiles:
            if ">>" in reaction_smiles:
                parts = reaction_smiles.split(">>")
                if len(parts) == 2:
                    reactants, products = parts
                    r_mol = parse_smiles(reactants)
                    p_mol = parse_smiles(products)
                    if r_mol and p_mol:
                        try:
                            detection = detect_reaction_types(reaction_smiles)
                            if detection.matches:
                                best = detection.matches[0]
                                self.validation_label.setText(f"OK Reaction: {best.name} detected")
                            else:
                                self.validation_label.setText("OK Reaction SMILES detected")
                        except Exception:
                            self.validation_label.setText("OK Reaction SMILES detected")
                        self.validation_label.setStyleSheet("color: #4cd137;")
                        return
            elif parse_smiles(reaction_smiles):
                self.validation_label.setText("OK Molecule SMILES detected")
                self.validation_label.setStyleSheet("color: #4cd137;")
                return

        if content.startswith("/"):
            self.validation_label.setText("")
            return

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

    def _extract_reaction_smiles(self, content: str) -> Optional[str]:
        text = (content or "").strip()
        if not text or ">>" not in text:
            return None

        if len(text) > 2 and text[0] == text[-1] and text[0] in {"`", '"', "'"}:
            text = text[1:-1].strip()

        match = REACTION_SMILES.search(text)
        if match:
            return match.group(1).strip()

        if " " not in text and "\t" not in text and "\n" not in text:
            return text
        return None

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
                self.append_chat("System", "Usage: /term <motif_id> <SMILES>")
                return
            motif_id = args[0]
            smiles = " ".join(args[1:]).strip()
        else:
            if not args:
                self.append_chat("System", "Usage: /terms <SMILES>")
                return
            motif_id = None
            smiles = " ".join(args).strip()

        try:
            from chemtools.featurizers.molecule import featurize_molecule

            payload = featurize_molecule(smiles)
            meta = payload.get("meta", {}) if isinstance(payload, dict) else {}
            error = meta.get("error") if isinstance(meta, dict) else None
            if error:
                self.append_chat("System", f"Motif analysis failed: {error}")
                return

            analyses = payload.get("analyses") or []
            if motif_id:
                analyses = [entry for entry in analyses if entry.get("compound_id") == motif_id]
                if not analyses:
                    self.append_chat("System", f"No motif matches for {motif_id}.")
                    return

            if not analyses:
                self.append_chat("System", "No motifs detected.")
                return

            lines = []
            for entry in analyses:
                compound_id = entry.get("compound_id", "unknown")
                groups = entry.get("nearby_groups") or []
                group_text = ", ".join(groups) if groups else "none"
                lines.append(f"- {compound_id}: {group_text}")

            header = "Nearby Groups (per motif):"
            self.append_chat("System", header + "\n" + "\n".join(lines))
        except Exception as exc:  # pragma: no cover - UI feedback
            self.append_chat("System", f"Term evaluation failed: {exc}")

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

        self._auto_render_from_response(response)

        self.append_chat("Agent", clean_response)
        self.append_log("Agent response received.")

    def _auto_render_from_response(self, text: str) -> None:
        if not text:
            return

        reaction_patterns = [
            r"(?:reaction|rxn)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[\S]+)",
            r"(?:SMILES|smiles)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[\S]+)",
            r"`([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)`",
            r"\"([A-Za-z0-9@+\-\[\]\(\)=#$\.>]+>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)\"",
            r"\b([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{10,}>>[A-Za-z0-9@+\-\[\]\(\)=#$\.>]+)\b",
        ]

        molecule_patterns = [
            r"(?:compound|molecule|structure)[\s:=]+([A-Za-z0-9@+\-\[\]\(\)=#$\.]{3,})",
            r"(?:SMILES|smiles)[\s:=]+([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{2,})",
            r"`([A-Za-z][A-Za-z0-9@+\-\[\]\(\)=#$\.]{5,})`",
        ]

        for pattern in reaction_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                smiles = match.group(1).strip()
                if ">>" in smiles and len(smiles) > 10:
                    self.render_image_from_smiles(
                        "reaction",
                        smiles,
                        source="auto-detected",
                        silent=True,
                    )
                    return

        for pattern in molecule_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                smiles = match.group(1).strip()
                if ">>" not in smiles and len(smiles) >= 5:
                    if (smiles[0].isalpha() or smiles[0] in "[(") and \
                       any(c in smiles for c in "()[]=#123456"):
                        self.render_image_from_smiles(
                            "molecule",
                            smiles,
                            source="auto-detected",
                            silent=True,
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

        self.agent = None
        self.append_log(
            "LLM mode set: "
            f"{self.llm_provider}/{self.llm_model} (T={self.llm_temperature:g})"
        )

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
                "- /taxonomy (show taxonomy status)\n"
                "- /terms <SMILES> (nearby groups per motif)\n"
                "- /term <motif_id> <SMILES> (nearby groups for one motif)\n"
                "- /term <motif_id> <SMILES> (nearby groups for one motif)\n"
                "- image <reaction|molecule> <SMILES>",
            )
            return False
        return True

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
            self.append_log("Dropped content loaded into input.")
            event.acceptProposedAction()
