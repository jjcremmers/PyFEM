# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from numpy import dot, zeros, array, cross
from math import sqrt

class DistributedLoad( Element ):
  
  def __init__ ( self, elnodes , props ):
    
    Element.__init__( self, elnodes , props )

    if self.rank == 2:
      self.dofTypes = [ 'u' , 'v' ]
    elif self.rank == 3:
      self.dofTypes = [ 'u' , 'v' , 'w' ]
      
    if hasattr(self,"trac"):
      self.trac = array(self.trac)
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def getExternalForce( self, elemdat ):
       
    sData = self.getShapeData( elemdat )
                       
    for iData in sData:
      N = self.getNmatrix( iData.h )
      
      trac = self.getTraction(iData.normal)
            
      elemdat.fint += dot(trac,N)*iData.weight      
      
    elemdat.fint *= self.loadFactor()
         
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( self.rank , self.rank*len(h) ) )

    for i,a in enumerate( h ):
      for j in list(range(self.rank)):
        N[j,self.rank*i+j] = a
    
    return N
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def getShapeData( self , elemdat ):
 
   crd  = elemdat.coords  
   nNod = crd.shape[0]
     
   if self.rank == 2:
     b = zeros(3)
     b[2] = 1.0
     
     if nNod == 2:
       sData = getElemShapeData( elemdat.coords , elemType = "Line2" )
       a     = self.getDirection(crd,1,0)
     elif nNod == 3:
       sData = getElemShapeData( elemdat.coords , elemType = "Line3" )
       a     = self.getDirection(crd,2,0)
     else:
       raise RuntimeError("The rank is 2, the number of nodes must be 2 or 3.")
   elif self.rank == 3:
     if nNod == 3:
       a     = self.getDirection(crd,1,0)
       b     = self.getDirection(crd,2,0)
     elif nNod == 4:
       sData = getElemShapeData( elemdat.coords , elemType = "Quad4" )  
       a     = self.getDirection(crd,1,0)
       b     = self.getDirection(crd,2,0)
     elif nNod == 6:
       sData = getElemShapeData( elemdat.coords , elemType = "Tria6" ) 
       a     = self.getDirection(crd,1,0)
       b     = self.getDirection(crd,2,0) 
     elif nNod == 8:
       sData = getElemShapeData( elemdat.coords , elemType = "Quad8" )  
       a     = self.getDirection(crd,1,0)
       b     = self.getDirection(crd,2,0)       
     else:
       raise RuntimeError("The rank is 3, the number of nodes must be 3, 4, 6 or 8.")
   else:
     raise RuntimeError("The element must be rank 3.")
     
   for iData in sData:
     iData.normal = cross(a,b)
     iData.normal *= 1.0/sqrt(dot(iData.normal,iData.normal))
     
   return sData     
   
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTraction( self , normal ):
  
    if hasattr(self,"pressure"):
      return self.solverStat.lam*normal * self.pressure    
    elif hasattr(self,"trac"):
      return self.solverStat.lam*self.trac
    else:
      raise RuntimeError("Define either pressure or trac")  
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getDirection( self , crd , i , j ):
    
    direc = zeros(3)
    
    rank = crd.shape[1]
    
    direc[:rank] = crd[i,:] - crd[j,:]
    
    return direc
    
       
              










           
