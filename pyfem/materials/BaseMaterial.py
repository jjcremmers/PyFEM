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

class BaseMaterial:

  def __init__ ( self, props ):

    for name,val in props:
      setattr( self, name, val )

    self.initHistory={}
    self.current = []
    self.iIter = -1

  def setIter( self, iIter ):
    
    self.iIter = iIter

  def setHistoryParameter( self , name , val ):

    if self.iIter == -1:
      self.initHistory[name]=val
      return
    
    if len(self.current) == self.iIter:
      self.current.append(self.initHistory.copy())
 
    self.current[self.iIter][name] = val
    
    return 
   
  def getHistoryParameter( self , name ):

    if len(self.history) == 0:
      return self.initHistory[name]
    else:
      return self.history[self.iIter][name]
    
  def commitHistory( self ):

    self.history = []

    for h in self.current:
      self.history.append(h)
  
