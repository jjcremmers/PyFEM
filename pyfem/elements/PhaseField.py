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
from pyfem.util.kinematics      import Kinematics
from numpy import zeros, dot, outer, ones , eye, ix_, linalg, tensordot 

import sys

class PhaseField( Element ):
  
  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.rank = props.rank
    self.k    = 1.0e-6

    if self.rank == 2:
      self.dofTypes = [ 'u' , 'v' , 'phase' ]
      self.nstr = 3
    elif self.rank == 3:
      self.dofTypes = [ 'u' , 'v' , 'w' , 'phase' ]
      self.nstr = 6

    self.kin = Kinematics(self.rank,self.nstr)
    
    self.hisOld = zeros(8)
    self.hisNew = zeros(8)

  def __type__ ( self ):
    return name
    
#-------------------------------------------------------------------------------
#  
#-------------------------------------------------------------------------------
    
  def strain2matrix ( self, strain ):
  
    '''Gives the strain in matrix format.'''    
    
    strainM = zeros( shape=( self.rank , self.rank ) )

    if self.rank == 2:
      strainM[0,0] = strain[0]
      strainM[1,1] = strain[1]
      strainM[0,1] = strain[2]
      strainM[1,0] = strain[2]

    elif self.rank == 3:
      strainM[0,0] = strain[0]
      strainM[1,1] = strain[1]
      strainM[2,2] = strain[2]
      strainM[1,2] = strain[3]
      strainM[0,2] = strain[4]
      strainM[0,1] = strain[5]
      strainM[2,1] = strain[3]
      strainM[2,0] = strain[4]
      strainM[1,0] = strain[5]
       
    return strainM
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getDecompEnergy ( self, strain ):
  
    '''Decomposes the Energy in a positive part due to tension and a negative part due to compression.'''    
  
    prinVal, prinVec = linalg.eig(strain)
    strainPos = zeros(shape=(self.rank,self.rank))
    strainNeg = zeros(shape=(self.rank,self.rank))
    
    for i in range(self.rank):
        strainPos += 0.5*(prinVal[i]+abs(prinVal[i]))*tensordot(prinVec[:,i],prinVec[:,i],0)
        strainNeg += 0.5*(prinVal[i]-abs(prinVal[i]))*tensordot(prinVec[:,i],prinVec[:,i],0)
      
    mu = self.matProps.E/(2*(1+self.matProps.nu))
    lame = (self.matProps.E*self.matProps.nu)/((1+self.matProps.nu)*(1-2*self.matProps.nu))
    
    energyPos = 0.5*lame*(0.5*(strain.trace()+abs(strain.trace())))**2 + mu*(strainPos*strainPos).trace()
    energyNeg = 0.5*lame*(0.5*(strain.trace()-abs(strain.trace())))**2 + mu*(strainNeg*strainNeg).trace()

    return energyPos, energyNeg

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    '''Calculates the tangent stiffness matrix and the internal force vector
       of the PhaseField model.'''
       
    sData = getElemShapeData( elemdat.coords )
    
    uDofs,pDofs = self.splitDofIDs( len(elemdat.coords) )
                 
    for iInt,iData in enumerate(sData):
      
      B = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( B , elemdat.state [uDofs] )
      self.kin.dstrain = dot ( B , elemdat.Dstate[uDofs] )
      
      phase     = dot( iData.h    , elemdat.state[pDofs] )
      gradPhase = dot( iData.dhdx.transpose() , elemdat.state[pDofs] )
      
      sigma,tang = self.mat.getStress( self.kin )

      energyPos, energyNeg = self.getDecompEnergy(self.strain2matrix(self.kin.strain))
      energy = 0.5*sum(self.kin.strain*sigma)

      if energyPos > self.hisOld[iInt]:
        self.hisNew[iInt] = energyPos
      else:
        self.hisNew[iInt] = self.hisOld[iInt]
            
      factor = (1.0-phase)**2+self.k
                   
      # -- Displacement contributions
      
      uStiff = dot ( B.transpose() , dot ( factor*tang , B ) ) 
      elemdat.stiff[ix_(uDofs,uDofs)] += uStiff * iData.weight
  
      dispForce = dot ( B.transpose() , factor*sigma )
      elemdat.fint[uDofs] += dispForce * iData.weight
      
      pStiff = (self.Gc/self.l0+2.0*self.hisNew[iInt])*outer(iData.h , iData.h )
      pStiff += self.Gc*self.l0*dot( iData.dhdx,iData.dhdx.transpose() )
      pStiff = iData.weight * pStiff
      
      elemdat.stiff[ix_(pDofs,pDofs)] += pStiff
      
      # -- Phase field contributions
      
      pfint =  self.Gc*self.l0*dot( iData.dhdx , gradPhase );
      pfint += self.Gc/self.l0*iData.h*phase;
      pfint += 2.0*( phase-1.0 ) * iData.h * self.hisNew[iInt]
            
      elemdat.fint[pDofs] += pfint * iData.weight
      
      # -- Coupling terms
      
      vecu = -2.0 * ( 1.0 - phase ) * dot( B.transpose() , sigma ) * iData.weight
      elemdat.stiff[ix_(uDofs,pDofs)] += outer( vecu , iData.h )       
      elemdat.stiff[ix_(pDofs,uDofs)] += outer( iData.h , vecu ) 
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
     
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
    
    '''Returns the internal force vector of the PhaseField model.'''
    
    sData = getElemShapeData( elemdat.coords )

    uDofs,pDofs = self.splitDofIDs( len(elemdat.coords) )

    for iInt,iData in enumerate(sData):
    
      B = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( B , elemdat.state [uDofs] )
      self.kin.dstrain = dot ( B , elemdat.Dstate[uDofs] )
      
      phase     = dot( iData.h    , elemdat.state[pDofs] )
      gradPhase = dot( iData.dhdx.transpose() , elemdat.state[pDofs] )
      
      sigma,tang = self.mat.getStress( self.kin )

      energyPos, energyNeg = self.getDecompEnergy(self.strain2matrix(self.kin.strain))
      energy = 0.5*sum(self.kin.strain*sigma)

      if energyPos > self.hisOld[iInt]:
        self.hisNew[iInt] = energyPos
      else:
        self.hisNew[iInt] = self.hisOld[iInt]
            
      factor = (1.0-phase)**2+self.k
                   
      # -- Displacement contributions
        
      dispForce = dot ( B.transpose() , factor*sigma )
      elemdat.fint[uDofs] += dispForce * iData.weight
           
      # -- Phase field contributions
      
      pfint =  self.Gc*self.l0*dot( iData.dhdx , gradPhase );
      pfint += self.Gc/self.l0*iData.h*phase;
      pfint += 2.0*( phase-1.0 ) * iData.h * self.hisNew[iInt]
            
      elemdat.fint[pDofs] += pfint * iData.weight
            
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def commit ( self, elemdat ):
  
    '''Copies the new history parameters (maximum internal energy for each 
       integration point) to the old history parameters.'''

    self.hisOld = self.hisNew
   
#-------------------------------------------------------------------------------
#  getBmatrix
#-------------------------------------------------------------------------------

  def getBmatrix( self , dphi ):

    '''Calculates the B-matrix (eps = B * u) vfor the mechanical part of the
       PhaseField model. The dimensions  of the B-matrix are determined by the
       rank of the problem and the length of the matrix containing derivatives
       of the shape functions (dphi)'''
       
    b = zeros( shape=( self.nstr , self.rank*len(dphi) ) )

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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def splitDofIDs( self , n ):
  
    '''Routine to split the dof IDs in two groups, one for the displacement 
       degrees of freedom, the second for the phase field degres of freedom. 
       n is the number of degrees of freedom in this model'''
    
    if self.rank == 2:
      if n == 3:
        return [0,1,3,4,6,7],[2,5,8]
      elif n == 4:
        return [0,1,3,4,6,7,9,10],[2,5,8,11]
    elif self.rank == 3:
      if n == 4:
        return [0,1,2,4,5,6,8,9,10,12,13,14],[3,7,11,15]
      elif n == 6:
        return [0,1,2,4,5,6,8,9,10,12,13,14,16,17,18,20,21,22],[3,7,11,15,19,23]
      elif n == 8:
        return [0,1,2,4,5,6,8,9,10,12,13,14,16,17,18,20,21,22,24,25,26,28,29,30],[3,7,11,15,19,23,27,31]
    else:
      print("Error")
