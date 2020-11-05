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

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from numpy import dot, outer, ix_

class ThermoSurface( Element ):
  
  def __init__ ( self, elnodes , props ):
  
    self.emissivity = 0.0
    self.convection = 0.0
    self.extTemp    = 0.0
    
    Element.__init__( self, elnodes , props )

    self.dofTypes = [ 'temp' ]
    self.Boltzman = 5.670373e-8
    self.extTemp4 = self.extTemp**4
        
  def __type__ ( self ):
    return name

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):
       
    sData = getElemShapeData( elemdat.coords , elemType = "Line2" )
                                  
    for iData in sData:
      temp     = sum( iData.h * elemdat.state )
                    
      elemdat.stiff += outer ( iData.h , iData.h ) * \
        ( self.convection + 4.0 * self.Boltzman * self.emissivity * temp**3 ) * iData.weight
           
      elemdat.fint += iData.h * ( self.convection * ( temp - self.extTemp ) + \
        self.Boltzman * self.emissivity * ( temp**4 - self.extTemp4 ) ) * iData.weight
           
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
     
    sData = getElemShapeData( elemdat.coords , elemType = "Line2" )
                       
    for iData in sData:         
      temp     = sum( iData.h * elemdat.state )
                               
      elemdat.fint += iData.h * ( self.convection * ( temp - self.extTemp ) + \
        self.Boltzman * self.emissivity * ( temp**4 - self.extTemp4 ) ) * iData.weight
