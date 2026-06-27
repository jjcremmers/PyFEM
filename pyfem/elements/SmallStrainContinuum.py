# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element
from pyfem.util.shapeFunctions  import getElemShapeData
from pyfem.util.kinematics      import Kinematics

import numpy as np

class SmallStrainContinuum( Element ):
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.rank = props.rank

    if self.rank == 2:
      self.dofTypes = [ 'u' , 'v' ]
      self.nstr = 3
      #self.outputLabels = ["s11","s22","s12"]
    elif self.rank == 3:
      self.dofTypes = [ 'u' , 'v' , 'w' ]
      self.nstr = 6
      #self.outputLabels = ["s11","s22","s33","s23","s13","s12"]

    self.kin = Kinematics(self.rank,self.nstr)

    #def __type__ ( self ):
    #return name

#------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    sData = getElemShapeData( elemdat.coords )

    #elemdat.outlabel.append(self.outputLabels)
    #elemdat.outdata  = zeros( shape=(len(elemdat.nodes),self.nstr) )

    for iData in sData:
      
      b = self.getBmatrix( iData.dhdx )

      self.kin.strain  = b @ elemdat.state
      self.kin.dstrain = b @ elemdat.Dstate
      
      sigma,tang = self.mat.getStress( self.kin )

      elemdat.stiff += b.T @ ( tang @ b ) * iData.weight
      elemdat.fint  += b.T @ sigma * iData.weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
      
      #import numpy as np
      #self.appendElementOutput( ["A","B"] , (self.iElm+1)*np.array([6,7]) )
     
#-------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    #elemdat.outlabel.append(self.outputLabels)
    #elemdat.outdata  = zeros( shape=(len(elemdat.nodes),self.nstr) )

    for iData in sData:
      b = self.getBmatrix( iData.dhdx )

      self.kin.strain  = b @ elemdat.state
      self.kin.dstrain = b @ elemdat.Dstate

      sigma,tang = self.mat.getStress( self.kin )

      elemdat.fint    += b.T @ sigma * iData.weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
 
#-------------------------------------------------------------------------------
    
  def getDissipation ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    for iData in sData:
      b = self.getBmatrix( iData.dhdx )

      self.kin.strain  = b @ elemdat.state
      self.kin.dstrain = b @ elemdat.Dstate

      self.mat.getStress( self.kin )

      self.kin.dgdstrain = np.zeros( 3 )
      self.kin.g = 0.0
            
      elemdat.fint += b.T @ self.kin.dgdstrain * iData.weight
      elemdat.diss += self.kin.g * iData.weight
           
#----------------------------------------------------------------------
    
  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    rho = elemdat.matprops.rho

    for iData in sData:
      N  = self.getNmatrix( iData.h )
      elemdat.mass += N.T @ N * rho * iData.weight
     
    elemdat.lumped = sum(elemdat.mass)
   
#--------------------------------------------------------------------------

  def getBmatrix( self , dphi ):

    b = np.zeros( shape=( self.nstr , self.dofCount() ) )

    if self.rank == 2:
      for i,dp in enumerate(dphi):
        b[0,i*2+0] = dp[0]
        b[1,i*2+1] = dp[1]
        b[2,i*2+0] = dp[1]
        b[2,i*2+1] = dp[0]
    elif self.rank == 3:
      for i,dp in enumerate(dphi):
        b[0,i*3+0] = dp[0]
        b[1,i*3+1] = dp[1]
        b[2,i*3+2] = dp[2]

        b[3,i*3+1] = dp[2]
        b[3,i*3+2] = dp[1]

        b[4,i*3+0] = dp[2]
        b[4,i*3+2] = dp[0]

        b[5,i*3+0] = dp[1]
        b[5,i*3+1] = dp[0]
   
    return b

#------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = np.zeros( shape=( self.rank , self.rank*len(h) ) )

    for i,a in enumerate( h ):
      for j in list(range(self.rank)):
        N[j,self.rank*i+j] = a
    
    return N
