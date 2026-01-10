# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from pyfem.util.kinematics      import Kinematics
from numpy import zeros, dot, outer, ones , eye
from math  import pi

class SmallStrainAxiSym( Element ):
  
  def __init__ ( self, elnodes , props ):
  
    Element.__init__( self, elnodes , props )

    if props.rank != 2:
      print("Error")
    
    self.dofTypes = [ 'u' , 'v' ]
   
    self.kin = Kinematics(3,6)

#------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    sData = getElemShapeData( elemdat.coords )

    for iData in sData:

      r      = dot( elemdat.coords[:,0] , iData.h )
      weight = 2.0*pi*r*iData.weight
             
      b = self.getBmatrix( iData.dhdx , iData.h , r )
      
      strain  = dot ( b , elemdat.state )
      #self.kin.dstrain = dot ( b , elemdat.Dstate )
      
      self.kin.strain[0] = strain[0]
      self.kin.strain[1] = strain[1]
      self.kin.strain[2] = strain[2]
      self.kin.strain[3] = 0.
      self.kin.strain[4] = 0.
      self.kin.strain[5] = strain[3]
      
      #self.kin.dstrain = dot ( b , elemdat.Dstate )
      
      sigma,tang = self.mat.getStress( self.kin )
      
      s4 = self.stress6to4( sigma )
      t4 = self.tang6to4  ( tang )

      elemdat.stiff += dot ( b.transpose() , dot ( t4 , b ) ) * weight
      elemdat.fint  += dot ( b.transpose() , s4 ) * weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
      
    
#-------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    for iData in sData:
    
      r      = dot( elemdat.coords[:,0] , iData.h )
      weight = 2.0*pi*r*iData.weight
      
      b = self.getBmatrix( iData.dhdx , iData.h , r )

      strain  = dot ( b , elemdat.state )
      #self.kin.dstrain = dot ( b , elemdat.Dstate )
      
      self.kin.strain[0] = strain[0]
      self.kin.strain[1] = strain[1]
      self.kin.strain[2] = strain[2]
      self.kin.strain[3] = 0.
      self.kin.strain[4] = 0.
      self.kin.strain[5] = strain[3]
      
      #self.kin.dstrain = dot ( b , elemdat.Dstate )

      sigma,tang = self.mat.getStress( self.kin )
      
      s4 = self.stress6to4( sigma )

      elemdat.fint    += dot ( b.transpose() , s4 ) * weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
      
#----------------------------------------------------------------------
    
  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    rho = elemdat.matprops.rho

    for iData in sData:
    
      r      = dot( elemdat.coords[0,:] , iData.h )
      weight = 2.0*pi*r*iData.weight
      
      N  = self.getNmatrix( iData.h )
      elemdat.mass += dot ( N.transpose() , N ) * rho * weight
     
    elemdat.lumped = sum(elemdat.mass)
   
#--------------------------------------------------------------------------

  def getBmatrix( self , dphi , phi , r ):

    b = zeros( shape=( 4 , self.dofCount() ) )

    invr = 1.0/r
    
    for i,(dp,p) in enumerate(zip(dphi,phi)):
      b[0,i*2+0] = dp[0]
      b[1,i*2+1] = dp[1]
      b[2,i*2+0] = p*invr
      b[3,i*2+0] = dp[1]
      b[3,i*2+1] = dp[0]
 
    return b

#------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( 2 , 2*len(h) ) )

    for i,a in enumerate( h ):
      for j in list(range(2)):
        N[j,2*i+j] = a
    
    return N
    
#

  def stress6to4( self , sigma ):
  
    s4 = zeros(4)
    
    s4[0] = sigma[0]
    s4[1] = sigma[1]
    s4[2] = sigma[2]
    s4[3] = sigma[5]
    
    return s4
    
  def tang6to4( self , tang ):
  
    t4 = zeros( shape=(4,4) )
    
    t4[:3,:3] = tang[:3,:3]
    t4[:3,3 ] = tang[:3,5]
    t4[3,:3]  = tang[5,:3]
    t4[3,3]   = tang[5,5]
    
    return t4
    
    
    
    
  
