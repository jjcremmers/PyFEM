# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

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
