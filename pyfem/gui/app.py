# SPDX-License-Identifier: MIT
# Copyright (c) 2011-2026 Joris J.C. Remmers

import re
import sys
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, Signal
from PySide6.QtGui import QAction, QFont, QIcon, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QStatusBar,
    QStyle,
    QTabWidget,
    QTextEdit,
    QToolBar,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pyfem.io.InputReader import InputRead
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver import Solver


class EmittingStream(QObject):
    text_written = Signal(str)

    def write(self, text):
        if text:
            self.text_written.emit(str(text))

    def flush(self):
        pass


class WorkerThread(QThread):
    finished_signal = Signal(bool, str)

    def __init__(self, input_file):
        super().__init__()
        self.input_file = input_file

    def run(self):
        try:
            props, globdat = InputRead(self.input_file)
            solver = Solver(props, globdat)
            output = OutputManager(props, globdat)

            while globdat.active and not self.isInterruptionRequested():
                solver.run(props, globdat)
                output.run(props, globdat)

            if self.isInterruptionRequested():
                globdat.active = False
                self.finished_signal.emit(False, "Analysis aborted by user.")
            else:
                self.finished_signal.emit(True, "Analysis completed.")

            globdat.close()
        except Exception as exc:  # pragma: no cover - GUI error path
            self.finished_signal.emit(False, f"Analysis failed: {exc}")


class MainWindow(QMainWindow):
    SECTION_PATTERN = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*\{")
    INPUT_PATTERN = re.compile(r'^\s*input\s*=\s*"([^"]+)"\s*;')
    OUTPUT_PATTERN = re.compile(r"outputModules\s*=\s*\[(.*?)\]\s*;", re.DOTALL)
    SOLVER_BLOCK_PATTERN = re.compile(r"solver\s*=\s*\{(.*?)\};", re.DOTALL)
    TYPE_PATTERN = re.compile(r'type\s*=\s*"([^"]+)"\s*;')

    def __init__(self):
        super().__init__()

        self.input_file = None
        self.worker = None
        self._stdout = sys.stdout
        self._stderr = sys.stderr

        self.setWindowTitle("PyFEM Workbench")
        self.resize(1680, 980)
        self.setMinimumSize(1280, 760)

        self.emitting_stream = EmittingStream()
        self.emitting_stream.text_written.connect(self.handle_output)
        sys.stdout = self.emitting_stream
        sys.stderr = self.emitting_stream

        self.create_actions()
        self.create_menu()
        self.create_toolbar()
        self.setStatusBar(QStatusBar(self))
        self.statusBar().showMessage("Ready")

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)

        self.ribbon = self.build_ribbon()
        layout.addWidget(self.ribbon)

        workspace = self.build_workspace()
        layout.addWidget(workspace, 1)

        self.apply_styles()
        self.populate_default_state()

    def closeEvent(self, event):  # pragma: no cover - Qt lifecycle
        sys.stdout = self._stdout
        sys.stderr = self._stderr
        super().closeEvent(event)

    def create_actions(self):
        self.new_action = QAction(
            self.style().standardIcon(QStyle.SP_FileIcon), "New", self
        )
        self.new_action.setShortcut(QKeySequence("Ctrl+N"))
        self.new_action.triggered.connect(self.new_model)

        self.open_action = QAction(
            self.style().standardIcon(QStyle.SP_DialogOpenButton), "Open", self
        )
        self.open_action.setShortcut(QKeySequence("Ctrl+O"))
        self.open_action.triggered.connect(self.load_file)

        self.save_action = QAction(
            self.style().standardIcon(QStyle.SP_DialogSaveButton), "Save", self
        )
        self.save_action.setShortcut(QKeySequence("Ctrl+S"))
        self.save_action.triggered.connect(self.save_file)

        self.run_action = QAction(
            self.style().standardIcon(QStyle.SP_MediaPlay), "Run", self
        )
        self.run_action.setShortcut(QKeySequence("Ctrl+R"))
        self.run_action.triggered.connect(self.execute_script)

        self.abort_action = QAction(
            self.style().standardIcon(QStyle.SP_BrowserStop), "Abort", self
        )
        self.abort_action.triggered.connect(self.abort_run)
        self.abort_action.setEnabled(False)

        self.about_action = QAction("About", self)
        self.about_action.triggered.connect(self.show_about)

        self.exit_action = QAction(QIcon.fromTheme("application-exit"), "Exit", self)
        self.exit_action.triggered.connect(self.close)

    def create_menu(self):
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("File")
        file_menu.addAction(self.new_action)
        file_menu.addAction(self.open_action)
        file_menu.addAction(self.save_action)
        file_menu.addSeparator()
        file_menu.addAction(self.exit_action)

        view_menu = menu_bar.addMenu("View")
        view_menu.addAction("Model Builder")
        view_menu.addAction("Settings")
        view_menu.addAction("Messages")

        run_menu = menu_bar.addMenu("Run")
        run_menu.addAction(self.run_action)
        run_menu.addAction(self.abort_action)

        help_menu = menu_bar.addMenu("Help")
        help_menu.addAction(self.about_action)

    def create_toolbar(self):
        toolbar = QToolBar("Quick Access", self)
        toolbar.setMovable(False)
        toolbar.setIconSize(toolbar.iconSize())
        toolbar.addAction(self.new_action)
        toolbar.addAction(self.open_action)
        toolbar.addAction(self.save_action)
        toolbar.addSeparator()
        toolbar.addAction(self.run_action)
        toolbar.addAction(self.abort_action)
        self.addToolBar(toolbar)

    def build_ribbon(self):
        ribbon = QTabWidget()
        ribbon.setDocumentMode(True)
        ribbon.setTabPosition(QTabWidget.North)
        ribbon.setMinimumHeight(142)

        ribbon.addTab(
            self.build_ribbon_page(
                [
                    ("Model Wizard", "Create a new model skeleton.", self.new_model),
                    ("Open", "Open a PyFEM project file.", self.load_file),
                    ("Save", "Write the current input preview.", self.save_file),
                    ("Run", "Start the active study.", self.execute_script),
                ]
            ),
            "Home",
        )
        ribbon.addTab(
            self.build_ribbon_page(
                [
                    ("Geometry", "Review the mesh input and topology file.", None),
                    ("Imports", "Track referenced data files.", None),
                    ("Selections", "Organize named model regions.", None),
                ]
            ),
            "Geometry",
        )
        ribbon.addTab(
            self.build_ribbon_page(
                [
                    ("Mesh", "Inspect the active discretization.", None),
                    ("Quality", "Check element and node counts.", None),
                    ("Preview", "Refresh the graphics dashboard.", self.refresh_dashboard),
                ]
            ),
            "Mesh",
        )
        ribbon.addTab(
            self.build_ribbon_page(
                [
                    ("Physics", "Browse element and material sections.", None),
                    ("Materials", "Inspect assigned material blocks.", None),
                    ("Loads", "Review constraints and load definitions.", None),
                ]
            ),
            "Physics",
        )
        ribbon.addTab(
            self.build_ribbon_page(
                [
                    ("Study", "Select solver and execution strategy.", None),
                    ("Parameters", "Prepare run-time parameters.", None),
                    ("Abort", "Request stop of a running analysis.", self.abort_run),
                ]
            ),
            "Study",
        )
        ribbon.addTab(
            self.build_ribbon_page(
                [
                    ("Results", "Review generated output modules.", None),
                    ("Messages", "Open the lower log panel.", self.focus_messages),
                    ("Report", "Summarize the loaded model.", self.refresh_dashboard),
                ]
            ),
            "Results",
        )

        return ribbon

    def build_ribbon_page(self, actions):
        page = QWidget()
        layout = QHBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        for label, description, callback in actions:
            button = QToolButton()
            button.setText(label)
            button.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
            button.setAutoRaise(False)
            button.setMinimumWidth(110)
            button.setMinimumHeight(88)
            button.setIcon(self.style().standardIcon(QStyle.SP_FileDialogDetailedView))
            button.setToolTip(description)
            if callback is not None:
                button.clicked.connect(callback)
            else:
                button.clicked.connect(
                    lambda _checked=False, text=label: self.statusBar().showMessage(
                        f"{text} panel is a design placeholder.", 3000
                    )
                )
            layout.addWidget(button)

        layout.addStretch(1)
        return page

    def build_workspace(self):
        splitter = QSplitter(Qt.Horizontal)

        left_panel = self.build_left_panel()
        center_panel = self.build_center_panel()
        right_panel = self.build_right_panel()

        splitter.addWidget(left_panel)
        splitter.addWidget(center_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([310, 930, 360])

        return splitter

    def build_left_panel(self):
        splitter = QSplitter(Qt.Vertical)

        builder_frame = self.panel_frame("Model Builder")
        builder_layout = builder_frame.layout()
        self.model_tree = QTreeWidget()
        self.model_tree.setHeaderHidden(True)
        self.model_tree.itemSelectionChanged.connect(self.update_selection_details)
        builder_layout.addWidget(self.model_tree)

        studies_frame = self.panel_frame("Study Queue")
        studies_layout = studies_frame.layout()
        self.study_list = QListWidget()
        self.study_list.addItems(
            [
                "Stationary / Nonlinear",
                "Mesh Preview",
                "Output Modules",
            ]
        )
        studies_layout.addWidget(self.study_list)

        splitter.addWidget(builder_frame)
        splitter.addWidget(studies_frame)
        splitter.setSizes([560, 180])
        return splitter

    def build_center_panel(self):
        splitter = QSplitter(Qt.Vertical)

        self.workspace_tabs = QTabWidget()
        self.workspace_tabs.setDocumentMode(True)

        graphics_page = QWidget()
        graphics_layout = QVBoxLayout(graphics_page)
        graphics_layout.setContentsMargins(0, 0, 0, 0)
        self.graphics_header = QLabel("Model Dashboard")
        self.graphics_header.setObjectName("workspaceHeader")
        graphics_layout.addWidget(self.graphics_header)

        self.graphics_view = QTextEdit()
        self.graphics_view.setReadOnly(True)
        self.graphics_view.setObjectName("graphicsView")
        graphics_layout.addWidget(self.graphics_view, 1)

        preview_page = QWidget()
        preview_layout = QVBoxLayout(preview_page)
        preview_layout.setContentsMargins(0, 0, 0, 0)
        preview_header = QLabel("Input Preview")
        preview_header.setObjectName("workspaceHeader")
        preview_layout.addWidget(preview_header)

        self.input_preview = QPlainTextEdit()
        editor_font = QFont("Courier New")
        editor_font.setStyleHint(QFont.Monospace)
        self.input_preview.setFont(editor_font)
        self.input_preview.setLineWrapMode(QPlainTextEdit.NoWrap)
        preview_layout.addWidget(self.input_preview, 1)

        self.workspace_tabs.addTab(graphics_page, "Graphics")
        self.workspace_tabs.addTab(preview_page, "Input File")

        messages_frame = self.panel_frame("Messages")
        messages_layout = messages_frame.layout()
        self.output_terminal = QTextEdit()
        self.output_terminal.setReadOnly(True)
        self.output_terminal.setMinimumHeight(200)
        messages_layout.addWidget(self.output_terminal)

        splitter.addWidget(self.workspace_tabs)
        splitter.addWidget(messages_frame)
        splitter.setSizes([660, 220])
        return splitter

    def build_right_panel(self):
        container = self.panel_frame("Settings")
        layout = container.layout()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(12)

        self.overview_group = QGroupBox("Selection")
        overview_layout = QFormLayout(self.overview_group)
        self.selection_name = QLineEdit()
        self.selection_name.setReadOnly(True)
        self.selection_type = QLineEdit()
        self.selection_type.setReadOnly(True)
        self.selection_source = QLineEdit()
        self.selection_source.setReadOnly(True)
        overview_layout.addRow("Node", self.selection_name)
        overview_layout.addRow("Role", self.selection_type)
        overview_layout.addRow("Source", self.selection_source)

        self.study_group = QGroupBox("Study Settings")
        study_layout = QFormLayout(self.study_group)
        self.study_solver = QLineEdit()
        self.study_solver.setReadOnly(True)
        self.study_mesh = QLineEdit()
        self.study_mesh.setReadOnly(True)
        self.study_outputs = QLineEdit()
        self.study_outputs.setReadOnly(True)
        study_layout.addRow("Solver", self.study_solver)
        study_layout.addRow("Mesh File", self.study_mesh)
        study_layout.addRow("Outputs", self.study_outputs)

        notes_group = QGroupBox("Actions")
        notes_layout = QVBoxLayout(notes_group)
        run_button = QPushButton("Run Active Study")
        run_button.clicked.connect(self.execute_script)
        save_button = QPushButton("Save Input Preview")
        save_button.clicked.connect(self.save_file)
        notes_layout.addWidget(run_button)
        notes_layout.addWidget(save_button)

        self.summary_card = QLabel()
        self.summary_card.setWordWrap(True)
        self.summary_card.setObjectName("summaryCard")
        self.summary_card.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)

        content_layout.addWidget(self.overview_group)
        content_layout.addWidget(self.study_group)
        content_layout.addWidget(notes_group)
        content_layout.addWidget(self.summary_card)
        content_layout.addStretch(1)

        scroll.setWidget(content)
        layout.addWidget(scroll)
        return container

    def panel_frame(self, title):
        frame = QFrame()
        frame.setObjectName("panelFrame")
        layout = QVBoxLayout(frame)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        label = QLabel(title)
        label.setObjectName("panelTitle")
        layout.addWidget(label)
        return frame

    def apply_styles(self):
        self.setStyleSheet(
            """
            QMainWindow {
                background: #dfe5eb;
            }
            QMenuBar, QToolBar, QTabWidget::pane {
                background: #eef2f6;
            }
            QTabBar::tab {
                background: #d7dfe8;
                color: #203243;
                padding: 8px 14px;
                border: 1px solid #b6c3d1;
                border-bottom: none;
                min-width: 88px;
            }
            QTabBar::tab:selected {
                background: #f8fafc;
            }
            #panelFrame {
                background: #f8fafc;
                border: 1px solid #bcc8d4;
                border-radius: 4px;
            }
            #panelTitle {
                color: #1e3247;
                font-size: 14px;
                font-weight: 700;
                letter-spacing: 0.4px;
                padding-bottom: 2px;
            }
            #workspaceHeader {
                color: #203243;
                font-size: 16px;
                font-weight: 700;
                padding: 4px 0 2px 2px;
            }
            #graphicsView {
                background: qlineargradient(
                    x1: 0, y1: 0, x2: 1, y2: 1,
                    stop: 0 #ffffff, stop: 1 #edf3f8
                );
                border: 1px solid #bcc8d4;
                color: #223648;
                font-size: 13px;
            }
            #summaryCard {
                background: #e9f0f6;
                border: 1px solid #bcc8d4;
                border-radius: 4px;
                color: #294054;
                padding: 10px;
            }
            QTreeWidget, QListWidget, QPlainTextEdit, QTextEdit, QLineEdit {
                background: #ffffff;
                border: 1px solid #bcc8d4;
                color: #1f3344;
                selection-background-color: #b9cde0;
                selection-color: #13232f;
            }
            QGroupBox {
                font-weight: 700;
                color: #213547;
                border: 1px solid #c7d1db;
                border-radius: 4px;
                margin-top: 10px;
                padding-top: 10px;
                background: #ffffff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 4px;
            }
            QPushButton, QToolButton {
                background: #e6edf4;
                border: 1px solid #b7c4d2;
                border-radius: 4px;
                padding: 6px 10px;
                color: #1d3143;
            }
            QPushButton:hover, QToolButton:hover {
                background: #d4e1ec;
            }
            """
        )

    def populate_default_state(self):
        self.graphics_view.setHtml(self.build_dashboard_html("No model loaded.", {}, []))
        self.input_preview.setPlainText(
            "# Open a .pro file to inspect and edit the PyFEM model input.\n"
        )
        self.selection_name.setText("Model Builder")
        self.selection_type.setText("Workbench")
        self.selection_source.setText("UI")
        self.study_solver.setText("Not loaded")
        self.study_mesh.setText("Not loaded")
        self.study_outputs.setText("Not loaded")
        self.summary_card.setText(
            "This GUI is structured as a simulation workbench: model tree on the "
            "left, active workspace in the center, and settings on the right."
        )
        self.populate_model_tree("PyFEM Model", {})

    def new_model(self):
        template = (
            'input = "mesh.dat";\n\n'
            "solver =\n{\n"
            '  type = "NonlinearSolver";\n'
            "  iterMax = 20;\n"
            "  tol = 1.0e-6;\n"
            "};\n\n"
            'outputModules = [ "vtk" , "output" ];\n'
        )
        self.input_file = None
        summary = self.parse_model_summary(template)
        self.populate_model_tree("Unsaved model", summary)
        self.update_study_fields("Unsaved model", summary)
        self.update_dashboard("Unsaved model", summary)
        self.workspace_tabs.setCurrentIndex(1)
        self.input_preview.setPlainText(template)
        self.statusBar().showMessage("New model template created.", 3000)
        self.output_terminal.append("Created new model template.")

    def load_file(self):
        input_file, _ = QFileDialog.getOpenFileName(
            self,
            caption="Open Input File",
            dir="examples",
            filter="PyFEM Input Files (*.pro);;All Files (*.*)",
        )

        if not input_file:
            return

        self.input_file = input_file
        try:
            text = Path(input_file).read_text(encoding="utf-8")
        except OSError as exc:
            QMessageBox.critical(self, "Open Failed", str(exc))
            return
        summary = self.parse_model_summary(text)

        self.input_preview.setPlainText(text)
        self.populate_model_tree(Path(input_file).name, summary)
        self.update_study_fields(input_file, summary)
        self.update_dashboard(input_file, summary)

        self.output_terminal.append(f"Loaded file: {input_file}")
        self.statusBar().showMessage(f"Loaded {Path(input_file).name}", 4000)

    def save_file(self):
        if self.input_file is None:
            target, _ = QFileDialog.getSaveFileName(
                self,
                caption="Save Input File",
                dir="examples",
                filter="PyFEM Input Files (*.pro)",
            )
            if not target:
                return
            self.input_file = target

        try:
            Path(self.input_file).write_text(
                self.input_preview.toPlainText(), encoding="utf-8"
            )
        except OSError as exc:
            QMessageBox.critical(self, "Save Failed", str(exc))
            return
        self.output_terminal.append(f"Saved file: {self.input_file}")
        self.statusBar().showMessage(f"Saved {Path(self.input_file).name}", 4000)

    def execute_script(self):
        if self.worker and self.worker.isRunning():
            self.output_terminal.append("An analysis is already running.")
            return

        if not self.input_file:
            self.output_terminal.append("Load or save an input file before running.")
            self.workspace_tabs.setCurrentIndex(1)
            return

        if Path(self.input_file).exists():
            self.save_file()

        self.worker = WorkerThread(self.input_file)
        self.worker.finished_signal.connect(self.analysis_finished)
        self.worker.start()

        self.run_action.setEnabled(False)
        self.abort_action.setEnabled(True)
        self.output_terminal.append(f"Running analysis: {self.input_file}")
        self.statusBar().showMessage("Analysis running...")

    def abort_run(self):
        if self.worker and self.worker.isRunning():
            self.worker.requestInterruption()
            self.output_terminal.append("Abort requested.")
            self.statusBar().showMessage("Abort requested...", 3000)

    def analysis_finished(self, success, message):
        self.run_action.setEnabled(True)
        self.abort_action.setEnabled(False)
        self.output_terminal.append(message)
        self.statusBar().showMessage(message, 5000)
        if success:
            self.workspace_tabs.setCurrentIndex(0)

    def handle_output(self, text):
        text = text.rstrip()
        if text:
            self.output_terminal.append(text)

    def focus_messages(self):
        self.output_terminal.setFocus()
        self.statusBar().showMessage("Messages panel focused.", 2000)

    def refresh_dashboard(self):
        text = self.input_preview.toPlainText()
        summary = self.parse_model_summary(text)
        source = self.input_file if self.input_file else "Unsaved model"
        self.update_dashboard(source, summary)
        self.update_study_fields(source, summary)
        self.statusBar().showMessage("Dashboard refreshed.", 2000)

    def populate_model_tree(self, root_name, summary):
        self.model_tree.clear()

        root = QTreeWidgetItem([root_name])
        self.model_tree.addTopLevelItem(root)

        groups = [
            ("Global Definitions", ["input", "parameters"]),
            ("Component 1", summary.get("sections", [])),
            ("Study", [summary.get("solver", "solver")]),
            ("Results", summary.get("outputs", []) or ["output"]),
        ]

        for label, items in groups:
            parent = QTreeWidgetItem([label])
            root.addChild(parent)
            for item in items:
                child = QTreeWidgetItem([item])
                parent.addChild(child)

        self.model_tree.expandAll()
        self.model_tree.setCurrentItem(root)

    def update_selection_details(self):
        item = self.model_tree.currentItem()
        if item is None:
            return

        name = item.text(0)
        parent = item.parent().text(0) if item.parent() else "Root"
        self.selection_name.setText(name)
        self.selection_type.setText(parent)
        self.selection_source.setText(Path(self.input_file).name if self.input_file else "Unsaved")

        self.summary_card.setText(
            f"Selected `{name}` from `{parent}`. Use the center workspace to inspect "
            "the input file and the lower messages pane to monitor analyses."
        )

    def update_study_fields(self, source, summary):
        mesh_file = summary.get("input", "Unknown")
        solver = summary.get("solver", "Not defined")
        outputs = ", ".join(summary.get("outputs", [])) or "Not defined"

        self.study_solver.setText(solver)
        self.study_mesh.setText(mesh_file)
        self.study_outputs.setText(outputs)

        self.selection_name.setText(Path(source).name if source else "Model Builder")
        self.selection_type.setText("PyFEM Project")
        self.selection_source.setText(str(source))

    def update_dashboard(self, source, summary):
        self.graphics_view.setHtml(self.build_dashboard_html(source, summary, summary.get("sections", [])))

    def build_dashboard_html(self, source, summary, sections):
        section_html = "".join(
            f"<li><b>{section}</b></li>" for section in sections[:10]
        ) or "<li><b>No sections detected</b></li>"
        outputs = ", ".join(summary.get("outputs", [])) or "No output modules detected"
        solver = summary.get("solver", "No solver block detected")
        mesh_file = summary.get("input", "No mesh/data file detected")

        return f"""
        <div style="font-family:Segoe UI,Arial,sans-serif; color:#223648;">
          <h2 style="margin-bottom:6px;">PyFEM Analysis Workspace</h2>
          <p style="margin-top:0;">
            COMSOL-inspired layout for editing, reviewing, and running a PyFEM model.
          </p>
          <table style="border-collapse:collapse; width:100%; margin:14px 0 18px 0;">
            <tr>
              <td style="padding:10px; border:1px solid #c5d1dc; background:#f5f8fb;"><b>Project</b><br>{source}</td>
              <td style="padding:10px; border:1px solid #c5d1dc; background:#f5f8fb;"><b>Solver</b><br>{solver}</td>
              <td style="padding:10px; border:1px solid #c5d1dc; background:#f5f8fb;"><b>Mesh/Data</b><br>{mesh_file}</td>
              <td style="padding:10px; border:1px solid #c5d1dc; background:#f5f8fb;"><b>Outputs</b><br>{outputs}</td>
            </tr>
          </table>
          <h3 style="margin-bottom:6px;">Model Sections</h3>
          <ul style="margin-top:0;">{section_html}</ul>
          <h3 style="margin-bottom:6px;">Workspace Intent</h3>
          <p style="margin-top:0;">
            Use the left tree to navigate the model structure, edit the source in the Input File tab,
            and run the active study while watching solver messages below.
          </p>
        </div>
        """

    def parse_model_summary(self, text):
        sections = []
        outputs = []
        solver = None
        input_file = None

        for line in text.splitlines():
            section_match = self.SECTION_PATTERN.match(line)
            if section_match:
                sections.append(section_match.group(1))

            input_match = self.INPUT_PATTERN.match(line)
            if input_match:
                input_file = input_match.group(1)

        solver_block = self.SOLVER_BLOCK_PATTERN.search(text)
        if solver_block:
            type_match = self.TYPE_PATTERN.search(solver_block.group(1))
            solver = type_match.group(1) if type_match else "solver"

        output_match = self.OUTPUT_PATTERN.search(text)
        if output_match:
            outputs = [
                token.strip().strip('"')
                for token in output_match.group(1).split(",")
            ]

        return {
            "input": input_file,
            "outputs": [item for item in outputs if item],
            "sections": sections,
            "solver": solver,
        }

    def show_about(self):
        QMessageBox.about(
            self,
            "About PyFEM Workbench",
            (
                "PyFEM Workbench\n\n"
                "A COMSOL-inspired desktop layout for browsing, editing, and "
                "running PyFEM analyses with PySide6."
            ),
        )


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
