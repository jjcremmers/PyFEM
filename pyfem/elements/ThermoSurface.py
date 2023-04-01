################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by                #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since        #
#  then augmented and maintained by Joris J.C. Remmers.                        #
#  All rights reserved.                                                        #
#                                                                              #
#  A github repository, with the most up to date version of the code,          #
#  can be found here:                                                          #
#     https://github.com/jjcremmers/PyFEM/                                     #
#     https://pyfem.readthedocs.io/                                            #	
#                                                                              #
#  The original code can be downloaded from the web-site:                      #
#     http://www.wiley.com/go/deborst                                          #
#                                                                              #
#  The code is open source and intended for educational and scientific         #
#  purposes only. If you use PyFEM in your research, the developers would      #
#  be grateful if you could cite the book.                                     #    
#                                                                              #
#  Disclaimer:                                                                 #
#  The authors reserve all rights but do not guarantee that the code is        #
#  free from errors. Furthermore, the authors shall not be liable in any       #
#  event caused by the use of the program.                                     #
################################################################################

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from numpy import dot, outer, ix_
from math import pi

class ThermoSurface( Element ):
  
  def __init__ ( self, elnodes , props ):
  
    self.emissivity   = 0.0
    self.convection   = 0.0
    self.extTemp      = 0.0
    self.axiSymmetric = False
    
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
           
    sData = self.getShapeData( elemdat )
                                      
    for iData in sData:
    
      if self.axiSymmetric:
        r      = dot( elemdat.coords[:,0] , iData.h )
        weight = 2.0*pi*r*iData.weight
      else:
        weight = iData.weight
      
      temp     = sum( iData.h * elemdat.state )
                    
      elemdat.stiff += outer ( iData.h , iData.h ) * \
        ( self.convection + 4.0 * self.Boltzman * self.emissivity * temp**3 ) * weight
           
      elemdat.fint += iData.h * ( self.convection * ( temp - self.extTemp ) + \
        self.Boltzman * self.emissivity * ( temp**4 - self.extTemp4 ) ) * weight
           
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
     
    sData = self.getShapeData( elemdat )
                       
    for iData in sData: 
    
      if self.axiSymmetric:
        r      = dot( elemdat.coords[:,0] , iData.h )
        weight = 2.0*pi*r*iData.weight
      else:
        weight = iData.weight        
        
      temp     = sum( iData.h * elemdat.state )
                               
      elemdat.fint += iData.h * ( self.convection * ( temp - self.extTemp ) + \
        self.Boltzman * self.emissivity * ( temp**4 - self.extTemp4 ) ) * weight
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def getShapeData( self , elemdat ):

   nNod = elemdat.coords.shape[0]
     
   if self.rank == 2:     
     if nNod == 2:
       return getElemShapeData( elemdat.coords , elemType = "Line2" )
     elif nNod == 3:
       return getElemShapeData( elemdat.coords , elemType = "Line3" )
     else:
       raise RuntimeError("The rank is 2, the number of nodes must be 2 or 3.")
   elif self.rank == 3:
     if nNod == 3:
       return getElemShapeData( elemdat.coords , elemType = "Tria3" )  
     elif nNod == 4:
       return getElemShapeData( elemdat.coords , elemType = "Quad4" )  
     elif nNod == 6:
       return getElemShapeData( elemdat.coords , elemType = "Tria6" ) 
     elif nNod == 8:
       return getElemShapeData( elemdat.coords , elemType = "Quad8" )  
     else:
       raise RuntimeError("The rank is 3, the number of nodes must be 3, 4, 6 or 8.")
   else:
     raise RuntimeError("The element must be rank 3.")
     
   return None   
