############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stable version can be downloaded from the web-site:          #
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/jjcremmers/PyFEM                                  #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #  
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################
import sys
from PyQt5 import QtGui,QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QAction, QFileDialog
from PyQt5.QtGui import QIcon

#from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
#import matplotlib.pyplot as plt
import os.path

from pyfem.io.InputReader   import InputRead
from pyfem.io.OutputManager import OutputManager

from pyfem.solvers.Solver   import Solver
 
class MyPopup(QWidget):

  def __init__(self):
    
    QWidget.__init__(self)

  def paintEvent(self, e):

    dc = QtGui.QPainter(self)
    dc.drawLine(0, 0, 100, 100)
    dc.drawLine(100, 0, 0, 100)
        
class XPyFEM(QMainWindow):
    
  def __init__(self):
    super(XPyFEM, self).__init__()
        
    self.initUI()
        
  def initUI(self):      

    #self.textEdit = QtGui.QTextEdit()
    #self.setCentralWidget(self.textEdit)
    self.statusBar()
    
    menubar = self.menuBar()
    fileMenu   = menubar.addMenu('File')
    editMenu   = menubar.addMenu('Edit')
    viewMenu   = menubar.addMenu('View')
    searchMenu = menubar.addMenu('Search')
    toolsMenu  = menubar.addMenu('Tools')
    helpMenu   = menubar.addMenu('Help')

    # Add exit button
    newFile = QAction(QIcon('new128.png'),'&New', self )
    newFile.setShortcut('Ctrl+N')
    newFile.setStatusTip('New project')
       
    openFile = QAction(QIcon('pyfem/qt/img/open128.png'),'&Open', self )
    openFile.setShortcut('Ctrl+O')
    openFile.setStatusTip('Open new File')
    openFile.triggered.connect(self.showDialog)
        
    saveFile = QAction(QIcon('pyfem/qt/img/save128.png'),'&Save', self )
    saveFile.setShortcut('Ctrl+S')
    saveFile.setStatusTip('Save File')
    saveFile.triggered.connect(self.showDialog)
        
    # Add exit button
    exitButton = QAction(QIcon('new128.png'),'&Exit', self)
    exitButton.setShortcut('Ctrl+Q')
    exitButton.setStatusTip('Exit application')
    exitButton.setIconText("Exit")
    exitButton.triggered.connect(self.close)
    fileMenu.addAction(exitButton)

    # Add exit button
    editProButton = QAction(QtGui.QIcon('pyfem/qt/img/edit128.png'),'&Edit', self)     
    editProButton.setStatusTip('EDit')
    editProButton.setIconText("Edit Pro")

    editProButton.triggered.connect(self.doit)
        
    # Add exit button
    runButton = QAction(QtGui.QIcon('pyfem/qt/img/run128.png'),'&Run', self)
    runButton.setShortcut('Ctrl+R')
    runButton.setStatusTip('Run simulation')
    runButton.setIconText("Run")
    runButton.triggered.connect(self.runIt)
        
    # Add exit button
    abortButton = QAction(QtGui.QIcon('pyfem/qt/img/abort128.png'),'&Exit', self)
    abortButton.setShortcut('Ctrl+C')
    abortButton.setStatusTip('Abort simulation')
    abortButton.setIconText("Abort")
    abortButton.triggered.connect(self.cancelIt)
 
    fileMenu.addAction(newFile)       
    fileMenu.addAction(openFile)
    fileMenu.addAction(saveFile)        
    fileMenu.addAction(exitButton)
      
    editMenu = menubar.addMenu('&Edit')
    editMenu.addAction(editProButton)       
#    editMenu.addAction(openFile)
#    editMenu.addAction(saveFile)        
#    editMenu.addAction(exitButton)
      
    runMenu = menubar.addMenu('&Run')
    runMenu.addAction(runButton)       
    runMenu.addAction(abortButton)
#    runMenu.addAction(saveFile)        
#    runMenu.addAction(exitButton)
        
    self.setGeometry(600, 400, 650, 400)
    self.setWindowTitle('PyFEM')

    self.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)        
    self.toolbar = self.addToolBar('Exit')
    self.toolbar.addAction(openFile)
    self.toolbar.addAction(runButton)
    self.toolbar.addAction(abortButton)

    """
    self.figure = plt.figure()
    self.canvas = FigureCanvas(self.figure)
    self.toolbar2 = NavigationToolbar(self.canvas, self)

    layout = QtGui.QVBoxLayout()
    layout.addWidget(self.toolbar)
    layout.addWidget(self.canvas)
    """
 
    self.setWindowIcon(QtGui.QIcon('pyfem/qt/img/pyfem_icon.png'))
    self.show()
               
  def showDialog(self):

    fname = QFileDialog.getOpenFileName(self, 'Open file','.pro')[0]

    os.chdir(os.path.dirname(str(fname)))
    self.props,self.globdat = InputRead( str(fname) )

    print("File ",os.path.basename(str(fname))," is ready. Press RUN to start the simulation.\n\n")

  def runIt(self):

    solver = Solver        ( self.props , self.globdat )
    output = OutputManager ( self.props , self.globdat )

    while self.globdat.active:
      solver.run( self.props , self.globdat )
      output.run( self.props , self.globdat )

    print("Simulation finished.")
    sys.exit()

  def cancelIt(self):

    print("Abort simulation.")

  def doit(self):
        
    print("Opening a new popup window...")
    self.w = MyPopup()
    self.w.setGeometry(100, 100, 400, 200)
    self.w.show()
        
def main():
    
  app = QApplication(sys.argv)
  ex  = XPyFEM()
  sys.exit(app.exec_())


if __name__ == '__main__':
  main()
