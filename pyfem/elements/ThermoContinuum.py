# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from numpy import zeros, dot, outer, eye

class ThermoContinuum( Element ):
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.rank = props.rank

    self.dofTypes = [ 'temp' ]
    
    self.D     = self.material.heatConductivity*eye(self.rank)
    self.capac = self.material.heatCapacity
    
    if self.rank == 2:    
      self.labels = [ "qx" , "qy" ]
    elif self.rank == 3:
      self.labels = [ "qx" , "qy" , "qz" ]
    
    self.transient = True
    self.theta = 1.0
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):
       
    sData = getElemShapeData( elemdat.coords )
    
    nDof = len(elemdat.coords)
    
    temp0 = elemdat.state - elemdat.Dstate
    
    if self.transient:
      ctt      = zeros(shape=(nDof,nDof))
      invdtime = 1.0/self.solverStat.dtime
                       
    for iInt,iData in enumerate(sData):        
      gradTemp = dot( iData.dhdx.transpose() , elemdat.state )
            
      elemdat.stiff += \
        dot ( iData.dhdx , dot( self.D , iData.dhdx.transpose() ) ) * iData.weight
     
      if self.transient:
        ctt += self.capac * outer( iData.h , iData.h ) * iData.weight
              
      self.appendNodalOutput( self.labels , dot(self.D,gradTemp) ) 
    
    if self.transient:  
      ktt0 = invdtime * ctt - elemdat.stiff * ( 1.0-self.theta )
      
      elemdat.stiff *= self.theta
      
      elemdat.stiff += invdtime * ctt 
        
    elemdat.fint += dot ( elemdat.stiff , elemdat.state )
      
    if self.transient:
      elemdat.fint += -dot ( ktt0 , temp0 )
     
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
     
    nDof = self.dofCount()
    elemdat.stiff = zeros( shape=(nDof,nDof) )
    
    self.getTangentStiffness( elemdat )
