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

import copy
from numpy import zeros

class BaseMaterial:

  def __init__ ( self, props ):

    self.numericalTangent = False
    self.storeOutputFlag = False
    
    for name,val in props:
      setattr( self, name, val )

    self.oldHistory = {}
    self.newHistory = {}

    self.outLabels  = []
    self.solverStat = props.solverStat
    
  def setHistoryParameter( self , name , val ):

    self.newHistory[name]=val
    return
       
  def getHistoryParameter( self , name ):

    if type(self.oldHistory[name]) == float:
      return self.oldHistory[name]
    else:
      return self.oldHistory[name].copy()
              
  def commitHistory( self ):

    self.oldHistory = copy.deepcopy(self.newHistory)
    
  def setOutputLabels ( self, labels ):
    
    self.outLabels = labels
    self.outData = zeros(len(self.outLabels))
    return
    
  def storeOutputs ( self, data ):
    if self.storeOutputFlag:
      self.outData = data
    return
  
