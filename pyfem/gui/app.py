# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import sys
#import subprocess
from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton, QTextEdit, QFileDialog, QVBoxLayout, QWidget, QToolBar, QMessageBox, QStyle, QSplitter, QHBoxLayout, QTabWidget, QLabel, QComboBox
from PySide6.QtCore import QProcess, Qt, QThread, Signal, QObject
from PySide6.QtGui import QIcon, QAction, QKeySequence

from pyfem.io.InputReader   import InputRead
from pyfem.io.OutputManager import OutputManager
from pyfem.solvers.Solver   import Solver

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class EmittingStream(QObject):
    text_written = Signal(str)

    def write(self, text):
        self.text_written.emit(str(text))

    def flush(self):
        pass 
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class WorkerThread(QThread):
    output_signal = Signal(str)

    def __init__(self, input_file):
        super().__init__()
        self.input_file = input_file

    def run(self):
        
        self.props,self.globdat = InputRead( self.input_file )
        
        self.solver = Solver        ( self.props , self.globdat )
        self.output = OutputManager ( self.props , self.globdat )

        while self.globdat.active:
            self.solver.run( self.props , self.globdat )
            self.output.run( self.props , self.globdat )

        self.globdat.close()        
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
                
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("PyFEM")
        
        self.create_menu()
        self.create_toolbar()

        '''
        self.layout = QVBoxLayout()
        '''        
        # Main container widget
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # Main horizontal layout with a splitter
        main_splitter = QSplitter(Qt.Horizontal)
        main_layout = QHBoxLayout(main_widget)
        main_layout.addWidget(main_splitter)

        # Left panel with tabs
        tab_widget = QTabWidget()
        tab_widget.setTabPosition(QTabWidget.West)
        tab_widget.setFixedWidth(400)

        # Adding tabs to the tab widget
        
        tabss = ["Mesh","Elements","Solver","Output"]
        
        for tabs in tabss:
            tab = QWidget()
            tab_layout = QVBoxLayout()
            tab.setLayout(tab_layout)            
            
            if tabs == "Solver":
                # Add a combo box for solver selection in the Solver tab
                solver_label = QLabel("Select a solver:")
                combo_box = QComboBox()
                combo_box.addItems(["SolverA", "SolverB", "SolverC"])
                
                # Add the label and combo box to the tab layout
                tab_layout.addWidget(solver_label)
                tab_layout.addWidget(combo_box)
                
                tab_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
                tab_layout.setSpacing(5)                
              
            else:
                # Default content for other tabs
                tab_label = QLabel(f"Content for Tab {tabs}")
                tab_layout.addWidget(tab_label)

            #tab.setLayout(tab_layout)
            tab_widget.addTab(tab, f"{tabs}")

        main_splitter.addWidget(tab_widget)

        # Right side layout with a vertical splitter
        right_splitter = QSplitter(Qt.Vertical)
        main_splitter.addWidget(right_splitter)

        # Image display area
        image_area = QLabel("Image Area")
        image_area.setMinimumHeight(300)
        image_area.setMinimumWidth(900)        
        image_area.setStyleSheet("background-color: white;")
        image_area.setAlignment(Qt.AlignCenter)
        right_splitter.addWidget(image_area)
             

        self.output_terminal = QTextEdit()
        self.output_terminal.setReadOnly(True)
        self.output_terminal.setMinimumHeight(200)        
        self.output_terminal.setStyleSheet("background-color: black; color: white; font-family: 'Courier New', monospace;")
        right_splitter.addWidget(self.output_terminal)

        '''
        container = QWidget()
        container.setLayout(self.layout)
        self.setCentralWidget(container)
        '''

        self.input_file = None
        
        # Set initial sizes of the split panels
        #main_splitter.setSizes([300, 900])
        #right_splitter.setSizes([600, 300])
        
        self.emitting_stream = EmittingStream()
        self.emitting_stream.text_written.connect(self.handle_output)
        sys.stdout = self.emitting_stream
        sys.stderr = self.emitting_stream        

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def create_menu(self):
        menu_bar = self.menuBar()
        
        file_menu  = menu_bar.addMenu("File")
        edit_menu  = menu_bar.addMenu("Edit")
        run_menu   = menu_bar.addMenu("Run")
        about_menu = menu_bar.addMenu("About")

        new_action = QAction(self.style().standardIcon(QStyle.SP_FileIcon), 'New', self)
        new_action.setShortcut(QKeySequence("Ctrl+N"))        
        new_action.triggered.connect(self.load_file)
        file_menu.addAction(new_action)

        open_action = QAction(self.style().standardIcon(QStyle.SP_DialogOpenButton), 'Open', self)
        open_action.setShortcut(QKeySequence("Ctrl+O"))        
        open_action.triggered.connect(self.load_file)
        file_menu.addAction(open_action)

        save_action = QAction(self.style().standardIcon(QStyle.SP_DialogSaveButton), 'Save', self)
        save_action.setShortcut(QKeySequence("Ctrl+S"))        
        save_action.triggered.connect(self.load_file)
        file_menu.addAction(save_action)
        
        exit_action = QAction(QIcon.fromTheme("application-exit"), "Exit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        execute_action = QAction(self.style().standardIcon(QStyle.SP_MediaPlay), "Execute", self)
        execute_action.setShortcut(QKeySequence("Ctrl+R"))  
        execute_action.triggered.connect(self.execute_script)
        run_menu.addAction(execute_action)

        abort_action = QAction(self.style().standardIcon(QStyle.SP_BrowserStop), "Abort", self)
        abort_action.setShortcut(QKeySequence("Ctrl+C"))          
        abort_action.triggered.connect(self.execute_script)
        run_menu.addAction(abort_action)
        
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about)
        about_menu.addAction(about_action)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def create_toolbar(self):
        toolbar = QToolBar("Main Toolbar")
        self.addToolBar(toolbar)

        new_action = QAction(self.style().standardIcon(QStyle.SP_FileIcon), 
                     'New', self)        
        new_action.triggered.connect(self.load_file)
        toolbar.addAction(new_action)

        open_action = QAction(self.style().standardIcon(QStyle.SP_DialogOpenButton),
                      'Open', self)        
        open_action.triggered.connect(self.load_file)
        toolbar.addAction(open_action)

        save_action = QAction(self.style().standardIcon(QStyle.SP_DialogSaveButton),
                      "Save", self)
        save_action.triggered.connect(self.execute_script)
        toolbar.addAction(save_action)        

        execute_action = QAction(self.style().standardIcon(QStyle.SP_MediaPlay),
                         "Execute", self)        
        execute_action.triggered.connect(self.execute_script)
        toolbar.addAction(execute_action)

        abort_action = QAction(self.style().standardIcon(QStyle.SP_BrowserStop),
                       "Abort", self)
        abort_action.triggered.connect(self.execute_script)
        toolbar.addAction(abort_action)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def load_file(self):
        file_dialog = QFileDialog()
        self.input_file, _ = file_dialog.getOpenFileName(self, 
            caption = "Open Input File", 
            dir     = "examples", 
            filter  = "PyFEM Input Files (*.pro);;All Files (*.*)")
        if self.input_file:
            self.output_terminal.append(f"Loaded file: {self.input_file}")

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def execute_script(self):
        if not self.input_file:
            self.output_terminal.append("Please load an input file first.")
            return
        
        self.worker = WorkerThread(self.input_file)
        self.worker.output_signal.connect(self.handle_output)
        self.worker.start()   
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
        
    def handle_output(self, text):
        self.output_terminal.append( text.rstrip() )
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
           
    def update_output(self):
        output = self.process.readAllStandardOutput().data().decode()
        error = self.process.readAllStandardError().data().decode()
        self.output_terminal.append(output)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def show_about(self):
        QMessageBox.about(self, "About", "PyFEM gui\n\nDeveloped with PySide6.")
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def main():
    app = QApplication(sys.argv)

    window = MainWindow()
    #window.resize(800, 600)
    window.show()

    sys.exit(app.exec())

