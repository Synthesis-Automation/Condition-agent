"""Main Qt window for the ChemTools assistant."""

from __future__ import annotations

import html
import traceback
from typing import List, Optional, Tuple

from PyQt6.QtCore import QObject, QRunnable, Qt, QThreadPool, pyqtSignal, QTimer
from PyQt6.QtWidgets import (
    QDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSplitter,
    QStatusBar,
    QTextEdit,
    QPlainTextEdit,
    QVBoxLayout,
    QWidget,
)
from PyQt6.QtGui import QFont

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
    RuleBuilderAutofillDialog,
    RuleBuilderDialog,
)


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


class ChemAssistantWindow(QMainWindow):
    """Main application window."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("ChemTools Assistant")
        self.resize(1200, 800)
        self._apply_default_font()

        self.agent = ChemToolsAgent(verbose=False)
        self.history: List[BaseMessage] = []
        self.constraint_spec = ConstraintSpec()
        self.constraint_text = ""

        self.thread_pool = QThreadPool()
        self._build_ui()
        self._update_constraint_summary()
        self.statusBar().showMessage("Ready")

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

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)

        # Chat log
        self.chat_view = QTextEdit()
        self.chat_view.setReadOnly(True)
        self.chat_view.setAcceptRichText(True)
        self.chat_view.setPlaceholderText("Conversation history will appear here.")
        self.chat_view.setMinimumHeight(480)
        splitter.addWidget(self.chat_view)

        # Info panel
        self.info_panel = QWidget()
        info_layout = QVBoxLayout(self.info_panel)

        self.constraint_label = QLabel()
        self.cache_label = QLabel()

        constraint_box = QGroupBox("Constraints")
        const_layout = QGridLayout(constraint_box)
        const_layout.addWidget(self.constraint_label, 0, 0, 1, 2)
        const_layout.addWidget(QLabel("Cache status:"), 1, 0)
        const_layout.addWidget(self.cache_label, 1, 1)
        info_layout.addWidget(constraint_box)

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setPlaceholderText("System log...")
        info_layout.addWidget(self.log_view)
        info_layout.addStretch()

        splitter.addWidget(self.info_panel)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        central_layout.addWidget(splitter)

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

        cache_clear_btn = QPushButton("Clear Cache")
        cache_clear_btn.clicked.connect(self.clear_cache)
        button_row.addWidget(cache_clear_btn)

        tools_btn = QPushButton("List Tools")
        tools_btn.clicked.connect(self.show_tool_summary)
        button_row.addWidget(tools_btn)

        builder_btn = QPushButton("Rule Builder Editor")
        builder_btn.clicked.connect(self.open_rule_builder_dialog)
        button_row.addWidget(builder_btn)

        autofill_btn = QPushButton("Autofill Draft")
        autofill_btn.clicked.connect(self.open_autofill_dialog)
        button_row.addWidget(autofill_btn)

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
        self.log_view.appendPlainText(text.strip())

    def send_message(self) -> None:
        content = self.input_edit.toPlainText().strip()
        if not content:
            return
        self.input_edit.clear()
        self.append_chat("You", content)
        self.statusBar().showMessage("Agent is thinking... |")
        self.set_input_enabled(False)
        self.start_spinner("Agent is thinking")

        worker = Worker(
            self.agent.chat,
            content,
            self.history,
            15,
            constraints=self.constraint_spec,
        )
        worker.signals.result.connect(self.handle_agent_result)
        worker.signals.error.connect(self.handle_agent_error)
        worker.signals.finished.connect(self.reset_after_task)
        self.thread_pool.start(worker)

    def handle_agent_result(self, payload: Tuple[str, List[BaseMessage]]) -> None:
        try:
            response, history = payload
        except ValueError:
            self.append_log("Unexpected agent payload.")
            return
        self.history = history
        self.append_chat("Agent", response)
        self.append_log("Agent response received.")

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
        if summary == "none":
            summary = "None"
        self.constraint_label.setText(summary or "None")

        prompt = format_constraints_for_prompt(self.constraint_spec)
        if prompt:
            self.append_log(f"Constraint prompt: {prompt}")

    def show_cache_stats(self) -> None:
        stats = recommendation_cache_stats()
        self.cache_label.setText(
            f"{stats['entries']} entries, {stats['hits']} hits"
        )
        QMessageBox.information(
            self,
            "Cache status",
            "\n".join(f"{k}: {v}" for k, v in stats.items()),
        )

    def clear_cache(self) -> None:
        clear_recommendation_cache()
        self.cache_label.setText("Cache cleared.")
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
