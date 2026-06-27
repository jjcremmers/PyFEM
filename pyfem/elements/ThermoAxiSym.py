# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from numpy import zeros, dot, outer, eye
from math import pi

class ThermoAxiSym( Element ):
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )
    
    self.rank = props.rank
        
    if self.rank != 2:   
      raise RuntimeError("This is an axisymmetric element. Please use an input mesh with rank 2.")    
               
    self.dofTypes = [ 'temp' ]
    
    self.D     = self.material.heatConductivity*eye(2)
    self.capac = self.material.heatCapacity
           
    self.labels = [ "qr" , "qz" ]
    
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
    
      r      = dot( elemdat.coords[:,0] , iData.h )
      weight = 2.0*pi*r*iData.weight
         
      gradTemp = dot( iData.dhdx.transpose() , elemdat.state )
            
      elemdat.stiff += \
        dot ( iData.dhdx , dot( self.D , iData.dhdx.transpose() ) ) * weight
     
      if self.transient:
        ctt += self.capac * outer( iData.h , iData.h ) * weight
              
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
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getBmatrix( self, dphi , phi , r ):
  
    b = zeros( shape=( 3 , self.dofCount() ) )

    invr = 1.0/r
    
    for i,(dp,p) in enumerate(zip(dphi,phi)):
      b[0,i] = dp[0]
      b[1,i] = dp[1]
      b[2,i] = 0 #p*invr
 
    return b
  
        
       

