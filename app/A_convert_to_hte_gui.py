import sys
import os
import json
import pandas as pd
from pathlib import Path
from typing import Optional, List

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from PyQt6 import QtWidgets, QtCore
from app.A_convert_to_hte_format import process_reaction_dataset

class ConversionWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal(bool, str)
    progress = QtCore.pyqtSignal(str)

    def __init__(
        self,
        input_path: str,
        output_path: str,
        drop_no_catalyst: bool,
        reagent_csv_path: str,
        new_reagents_path: str,
    ):
        super().__init__()
        self.input_path = input_path
        self.output_path = output_path
        self.drop_no_catalyst = drop_no_catalyst
        self.reagent_csv_path = reagent_csv_path
        self.new_reagents_path = new_reagents_path

    def run(self):
        try:
            # Custom stdout capturing to redirect prints to the GUI log
            class StdoutRedirector:
                def __init__(self, signal):
                    self.signal = signal
                def write(self, text):
                    if text.strip():
                        self.signal.emit(text.strip())
                def flush(self):
                    pass

            old_stdout = sys.stdout
            sys.stdout = StdoutRedirector(self.progress)
            
            try:
                process_reaction_dataset(
                    self.input_path, 
                    self.output_path, 
                    drop_no_catalyst=self.drop_no_catalyst,
                    reagent_csv_path=self.reagent_csv_path,
                    new_reagents_path=self.new_reagents_path,
                )
                self.finished.emit(True, f"Successfully processed {self.input_path}")
            finally:
                sys.stdout = old_stdout
                
        except Exception as e:
            self.finished.emit(False, f"Error: {str(e)}")

class HTEConverterWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HTE Dataset Converter (Qt6)")
        self.resize(700, 500)
        
        # UI Elements
        self.dataset_dropdown = QtWidgets.QComboBox()
        self.source_dir_edit = QtWidgets.QLineEdit()
        self.source_dir_edit.setPlaceholderText("Source folder containing *.csv files")
        self.input_edit = QtWidgets.QLineEdit()
        self.output_edit = QtWidgets.QLineEdit()
        self.drop_catalyst_check = QtWidgets.QCheckBox("Drop reactions without a catalyst")
        self.drop_catalyst_check.setChecked(True)
        
        self.btn_run = QtWidgets.QPushButton("Start Conversion")
        self.btn_quit = QtWidgets.QPushButton("Quit")
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        
        self._setup_layout()
        self._load_datasets()
        
        # Connections
        self.dataset_dropdown.currentIndexChanged.connect(self._on_dataset_changed)
        self.source_dir_edit.editingFinished.connect(self._load_datasets)
        self.btn_run.clicked.connect(self.run_processing)
        self.btn_quit.clicked.connect(self.close)
        
        # Threading state
        self.thread = None
        self.worker = None

    def _setup_layout(self):
        layout = QtWidgets.QVBoxLayout(self)
        
        title = QtWidgets.QLabel("HTE Dataset Converter")
        title.setStyleSheet("font-size: 18px; font-weight: bold; margin: 10px;")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)
        
        form = QtWidgets.QFormLayout()
        
        source_layout = QtWidgets.QHBoxLayout()
        source_layout.addWidget(self.source_dir_edit)
        source_browse_btn = QtWidgets.QPushButton("Browse")
        source_browse_btn.setFixedWidth(80)
        source_browse_btn.clicked.connect(self._choose_source_dir)
        source_layout.addWidget(source_browse_btn)
        form.addRow("Source Folder:", source_layout)

        # Dataset selector
        ds_layout = QtWidgets.QHBoxLayout()
        ds_layout.addWidget(self.dataset_dropdown)
        refresh_btn = QtWidgets.QPushButton("Refresh")
        refresh_btn.setFixedWidth(80)
        refresh_btn.clicked.connect(self._load_datasets)
        ds_layout.addWidget(refresh_btn)
        form.addRow("Available Datasets:", ds_layout)
        
        form.addRow("Input CSV:", self.input_edit)
        form.addRow("Output CSV:", self.output_edit)
        form.addRow("", QtWidgets.QLabel("Output includes spectator_groups column."))
        form.addRow("", self.drop_catalyst_check)
        
        layout.addLayout(form)
        
        layout.addWidget(QtWidgets.QLabel("Processing Log:"))
        layout.addWidget(self.log)
        layout.addWidget(self.progress_bar)
        
        btns = QtWidgets.QHBoxLayout()
        btns.addStretch()
        btns.addWidget(self.btn_run)
        btns.addWidget(self.btn_quit)
        layout.addLayout(btns)

    def _load_datasets(self):
        self.dataset_dropdown.clear()
        self.dataset_dropdown.addItem("-- Select or enter custom path --", None)

        source_text = self.source_dir_edit.text().strip()
        dataset_dir = Path(source_text) if source_text else PROJECT_ROOT / "data" / "reaction_dataset"
        self.source_dir_edit.setText(str(dataset_dir))

        if dataset_dir.exists():
            files = sorted(p for p in dataset_dir.rglob("*.csv") if p.is_file())
            for f in files:
                self.dataset_dropdown.addItem(f.name, str(f))

    def _on_dataset_changed(self, index):
        if index <= 0:
            return
            
        data_path = self.dataset_dropdown.currentData()
        if not data_path:
            return

        input_path = Path(data_path)
        output_path = PROJECT_ROOT / "data" / "HTE_db" / "datasets" / f"{input_path.stem}_canonical.csv"

        self.input_edit.setText(str(input_path))
        self.output_edit.setText(str(output_path))

    def _choose_source_dir(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select Source Folder",
            str(PROJECT_ROOT),
        )
        if path:
            self.source_dir_edit.setText(path)
            self._load_datasets()

    def log_msg(self, text: str):
        self.log.appendPlainText(text)

    def run_processing(self):
        input_path = self.input_edit.text().strip()
        output_path = self.output_edit.text().strip()
        
        if not input_path or not os.path.exists(input_path):
            QtWidgets.QMessageBox.critical(self, "Error", "Invalid input path.")
            return
        if not output_path:
            QtWidgets.QMessageBox.critical(self, "Error", "Output path is required.")
            return

        self.setEnabled(False)
        self.log.clear()
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0) # Indeterminate
        
        # Start worker thread
        self.thread = QtCore.QThread()
        self.worker = ConversionWorker(
            input_path, 
            output_path, 
            self.drop_catalyst_check.isChecked(),
            str(PROJECT_ROOT / "data" / "reagent_db" / "reagents.csv"),
            str(PROJECT_ROOT / "data" / "reagent_db" / "new_reagents.csv"),
        )
        self.worker.moveToThread(self.thread)
        
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self.log_msg)
        self.worker.finished.connect(self.on_finished)
        
        self.thread.start()

    def on_finished(self, success: bool, message: str):
        self.thread.quit()
        self.thread.wait()
        
        self.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        if success:
            QtWidgets.QMessageBox.information(self, "Success", message)
            self.log_msg("\n--- Processing Complete ---")
        else:
            QtWidgets.QMessageBox.warning(self, "Error", message)

def main():
    app = QtWidgets.QApplication(sys.argv)
    window = HTEConverterWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
