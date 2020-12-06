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

from pyfem.util.dataStructures import Properties

class MaterialManager ( list ):

  def __init__ ( self, matProps ):

    matType = matProps.type

    self.material = getattr(__import__('pyfem.materials.'+matType , globals(), locals(), matType , 0 ), matType )
    
    self.matlist     = []
    self.matProps    = matProps
    self.iSam        = -1
    self.failureFlag = False
    
    if hasattr(matProps,'failureType'):
    
      failureType = matProps.failureType
      
      failure = getattr(__import__('pyfem.materials.'+failureType , \
        globals(), locals(), failureType , 0 ), failureType )
 
      self.failure = failure( matProps )
      self.failureFlag = True

  def reset( self ):

    self.iSam  = -1

  def getStress ( self, kinematic , iSam = -1 ):

    if iSam == -1:
      self.iSam += 1
    else:
      self.iSam = iSam
            
    while self.iSam >= len(self.matlist):
      self.matlist.append(self.material( self.matProps ))
        
    self.mat = self.matlist[self.iSam]
     
    result = self.mat.getStress( kinematic )
    
    if self.failureFlag:
      self.failure.check(result[0],kinematic)
      
    return result
    
  def outLabels( self ):
    return self.mat.outLabels

  def outData( self ):
    return self.mat.outData

  def getHistory( self , label ):
    return self.mat.getHistoryParameter( label )

  def commitHistory( self ):
    for mat in self.matlist:
      mat.commitHistory()
