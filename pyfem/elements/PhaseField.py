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
from pyfem.util.kinematics      import Kinematics
from numpy import zeros, dot, outer, ones , eye, ix_ 

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
      print("Error")

    self.kin = Kinematics(self.rank,self.nstr)
    
    self.hisOld = zeros(4)
    self.hisNew = zeros(4)

  def __type__ ( self ):
    return name

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    print(self.rank)
    sData = getElemShapeData( elemdat.coords )
    
    uDofs,pDofs = self.splitDofIDs( len(elemdat.coords) )
                 
    for iInt,iData in enumerate(sData):
      
      B = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( B , elemdat.state [uDofs] )
      self.kin.dstrain = dot ( B , elemdat.Dstate[uDofs] )
      
      phase     = dot( iData.h    , elemdat.state[pDofs] )
      gradPhase = dot( iData.dhdx.transpose() , elemdat.state[pDofs] )
      
      sigma,tang = self.mat.getStress( self.kin )

      energy = 0.5*sum(self.kin.strain*sigma)

      if energy > self.hisOld[iInt]:
        self.hisNew[iInt] = energy
      else:
        self.hisNew[iInt] = self.hisOld[iInt]
            
      factor = 1.0-phase*phase+self.k
                   
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
      
      # Coupling terms need TOBEIMPLEMENTED
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
     
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    uDofs,pDofs = self.splitDofIDs( len(elemdat.coords) )

    for iInt,iData in enumerate(sData):
    
      B = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( B , elemdat.state [uDofs] )
      self.kin.dstrain = dot ( B , elemdat.Dstate[uDofs] )
      
      phase     = dot( iData.h    , elemdat.state[pDofs] )
      gradPhase = dot( iData.dhdx.transpose() , elemdat.state[pDofs] )
      
      sigma,tang = self.mat.getStress( self.kin )

      energy = 0.5*sum(self.kin.strain*sigma)

      if energy > self.hisOld[iInt]:
        self.hisNew[iInt] = energy
      else:
        self.hisNew[iInt] = self.hisOld[iInt]
            
      factor = 1.0-phase*phase+self.k
                   
      # -- Displacement contributions
        
      dispForce = dot ( B.transpose() , factor*sigma )
      elemdat.fint[uDofs] += dispForce * iData.weight
           
      # -- Phase field contributions
      
      pfint =  self.Gc*self.l0*dot( iData.dhdx , gradPhase );
      pfint += self.Gc/self.l0*iData.h*phase;
      pfint += 2.0*( phase-1.0 ) * iData.h * self.hisNew[iInt]
            
      elemdat.fint[pDofs] += pfint * iData.weight
      
      # -- Coupling terms
  
      # Coupling terms need TOBEIMPLEMENTED
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )

#-------------------------------------------------------------------------------
#  
#-------------------------------------------------------------------------------
    
  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    rho = elemdat.matprops.rho

    for iData in sData:
      N  = self.getNmatrix( iData.h )
      elemdat.mass += dot ( N.transpose() , N ) * rho * iData.weight
     
    elemdat.lumped = sum(elemdat.mass)
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def commit ( self, elemdat ):

    self.hisOld = self.hisNew
   
#-------------------------------------------------------------------------------
# Calculates the B matrix
#-------------------------------------------------------------------------------

  def getBmatrix( self , dphi ):

    b = zeros( shape=( self.nstr , self.rank*len(dphi) )

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

  def getNmatrix( self , h ):

    N = zeros( shape=( self.rank , self.rank*len(h) ) )

    for i,a in enumerate( h ):
      for j in list(range(self.rank)):
        N[j,self.rank*i+j] = a
    
    return N
    
#-------------------------------------------------------------------------------
#  Routine to split the dof IDs in two groups, one for the displacement degrees 
#    of freedom, the second for the phase field degres of freedom
#-------------------------------------------------------------------------------

  def splitDofIDs( self , n ):
  
    if self.rank == 2:
      if n == 3:
        return [0,1,3,4,6,7],[2,5,8]
      elif n == 4:
        return [0,1,3,4,6,7,9,10],[2,5,8,11]
    elif self.rank == 3:
      if n == 8:
        return [0,1,2,4,5,6,8,9,10,12,13,14,16,17,18,20,21,22,24,25,26,28,29,30],[3,7,11,15,19,23,27,31]
    else:
      print("Error")
