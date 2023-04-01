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

from numpy import zeros, dot, outer, ones, eye, sqrt, reshape

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class FiniteStrainContinuum( Element ):
  
  def __init__ ( self, elnodes , props ):
  
    self.method = "TL"
    
    Element.__init__( self, elnodes , props )

    self.rank = props.rank

    if self.rank == 2:
      self.dofTypes = [ 'u' , 'v' ]
      self.nstr = 3
    elif self.rank == 3:
      self.dofTypes = [ 'u' , 'v' , 'w' ]
      self.nstr = 6

    self.kin = Kinematics(self.rank,self.nstr)
    
  def __type__ ( self ):
    return name

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    if self.method == "TL":
      return self.getTLTangentStiffness( elemdat )
    elif self.method == "UL":
      return self.getULTangentStiffness( elemdat )
    else:
      raise RuntimeError("Please define total or updated Lagrange (method = 'TL' or 'UL'.")
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTLTangentStiffness ( self, elemdat ):

    sData = getElemShapeData( elemdat.coords )
   
    for iData in sData:

      self.kin = self.getKinematics( iData.dhdx , elemdat ) 
      B        = self.getBmatrix   ( iData.dhdx , self.kin.F )
      
      sigma,tang = self.mat.getStress( self.kin )
        
      elemdat.stiff += dot ( B.transpose() , dot ( tang , B ) ) * iData.weight

      T   = self.stress2matrix( sigma )
      Bnl = self.getBNLmatrix ( iData.dhdx )
   
      elemdat.stiff += dot ( Bnl.transpose() , dot( T , Bnl ) ) * iData.weight
      elemdat.fint  += dot ( B.transpose() , sigma ) * iData.weight
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getULTangentStiffness ( self, elemdat ):
  
    state0  = elemdat.state - elemdat.Dstate
    curCrds = elemdat.coords + reshape(state0,elemdat.coords.shape)
    
    sData0 = getElemShapeData( elemdat.coords )
    sDataC = getElemShapeData( curCrds )

       
    for iOrig,iCurr in zip(sData0,sDataC):

      self.kin = self.getKinematics( iOrig.dhdx , elemdat ) 
      B        = self.getULBmatrix ( iCurr.dhdx )
            
      sigma,tang = self.mat.getStress( self.kin )
        
      elemdat.stiff += dot ( B.transpose() , dot ( tang , B ) ) * iCurr.weight

      T   = self.stress2matrix( sigma )
      Bnl = self.getBNLmatrix ( iCurr.dhdx )
   
      elemdat.stiff += dot ( Bnl.transpose() , dot( T , Bnl ) ) * iCurr.weight
      elemdat.fint  += dot ( B.transpose() , sigma ) * iCurr.weight
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
   
    if self.method == "TL":
      return self.getTLInternalForce( elemdat )
    elif self.method == "UL":
      return self.getULInternalForce( elemdat )
    else:
      raise RuntimeError("Please define total or updated Lagrange (method = 'TL' or 'UL'.")
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
  def getTLInternalForce ( self, elemdat ):

    sData = getElemShapeData( elemdat.coords )
       
    for iData in sData:

      self.kin = self.getKinematics( iData.dhdx , elemdat ) 
      B        = self.getBmatrix   ( iData.dhdx , self.kin.F )
      
      sigma,tang = self.mat.getStress( self.kin )
        
      elemdat.fint  += dot ( B.transpose() , sigma ) * iData.weight
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getULInternalForce ( self, elemdat ):
  
    state0  = elemdat.state - elemdat.Dstate
    curCrds = elemdat.coords + reshape(state0,elemdat.coords.shape)
 
    sData0 = getElemShapeData( elemdat.coords )
    sDataC = getElemShapeData( curCrds )
    
    for iOrig,iCurr in zip(sData0,sDataC):

      self.kin = self.getKinematics( iOrig.dhdx , elemdat ) 
      B        = self.getULBmatrix ( iCurr.dhdx )
            
      sigma,tang = self.mat.getStress( self.kin )
             
      elemdat.fint  += dot ( B.transpose() , sigma ) * iCurr.weight
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
      
#-------------------------------------------------------------------------------
    
  def getDissipation ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    for iData in sData:
      self.kin = self.getKinematics( iData.dhdx , elemdat ) 
      B        = self.getBmatrix   ( iData.dhdx , self.kin.F )

      self.mat.getStress( self.kin )

      self.kin.dgdstrain = zeros( 3 )
      self.kin.g = 0.0
      
      elemdat.fint += dot ( B.transpose() , self.kin.dgdstrain ) * iData.weight
      elemdat.diss += self.kin.g * iData.weight 
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    rho = elemdat.matprops.rho

    for iData in sData:
      N  = self.getNmatrix( iData.h )
      elemdat.mass += dot ( N.transpose() , N ) * rho * iData.weight
     
    elemdat.lumped = sum(elemdat.mass)
            
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getKinematics( self , dphi , elemdat ):
  
    kin = Kinematics(self.rank,self.nstr)
    
    elstate  = elemdat.state
    elstate0 = elstate - elemdat.Dstate
    
    kin.F = eye(self.rank)
    kin.F0 = eye(self.rank)
  
    for i in range(len(dphi)):
      for j in range(self.rank):
        for k in range(self.rank):
          kin.F[j,k] += dphi[i,k]*elstate[self.rank*i+j]
          kin.F0[j,k] += dphi[i,k]*elstate0[self.rank*i+j]

    kin.E = 0.5*(dot(kin.F.transpose(),kin.F)-eye(self.rank))

    kin.strain[0] = kin.E[0,0]
    kin.strain[1] = kin.E[1,1]

    if self.rank == 2:
      kin.strain[2] = 2.0*kin.E[0,1]
    elif self.rank == 3:
      kin.strain[2] =     kin.E[2,2]
      kin.strain[3] = 2.0*kin.E[1,2]
      kin.strain[4] = 2.0*kin.E[0,2]
      kin.strain[5] = 2.0*kin.E[0,1]
    
    return kin

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBmatrix( self , dphi , F ):

    B = zeros( shape=(self.nstr, self.rank*len(dphi) ) )

    if self.rank == 2:
      for i,dp in enumerate( dphi ):
        B[0,2*i  ] = dp[0]*F[0,0]
        B[0,2*i+1] = dp[0]*F[1,0]
	
        B[1,2*i  ] = dp[1]*F[0,1]
        B[1,2*i+1] = dp[1]*F[1,1]

        B[2,2*i  ] = dp[1]*F[0,0]+dp[0]*F[0,1]
        B[2,2*i+1] = dp[0]*F[1,1]+dp[1]*F[1,0]
    elif self.rank == 3:
      for i,dp in enumerate( dphi ):
        B[0,3*i  ] = dp[0]*F[0,0]
        B[0,3*i+1] = dp[0]*F[1,0]
        B[0,3*i+2] = dp[0]*F[2,0]
	
        B[1,3*i  ] = dp[1]*F[0,1]
        B[1,3*i+1] = dp[1]*F[1,1]
        B[1,3*i+2] = dp[1]*F[2,1]

        B[2,3*i  ] = dp[2]*F[0,2]
        B[2,3*i+1] = dp[2]*F[1,2]
        B[2,3*i+2] = dp[2]*F[2,2]

        B[3,3*i  ] = dp[1]*F[0,2]+dp[2]*F[0,1]
        B[3,3*i+1] = dp[1]*F[1,2]+dp[2]*F[1,1]
        B[3,3*i+2] = dp[1]*F[2,2]+dp[2]*F[2,1]

        B[4,3*i  ] = dp[2]*F[0,0]+dp[0]*F[0,2]
        B[4,3*i+1] = dp[2]*F[1,0]+dp[0]*F[1,2]
        B[4,3*i+2] = dp[2]*F[2,0]+dp[0]*F[2,2]

        B[5,3*i  ] = dp[0]*F[0,1]+dp[1]*F[0,0]
        B[5,3*i+1] = dp[0]*F[1,1]+dp[1]*F[1,0]
        B[5,3*i+2] = dp[0]*F[2,1]+dp[1]*F[2,0]
 
    return B
    
 #------------------------------------------------------------------------------
 #
 #------------------------------------------------------------------------------
    
  def getULBmatrix( self , dphi  ):

    B = zeros( shape=(self.nstr, self.rank*len(dphi) ) )

    if self.rank == 2:
      for iNel,dp in enumerate( dphi ):
        i = 2 * iNel
        
        B[0,i  ] = dp[0]	
        B[1,i+1] = dp[1]

        B[2,i  ] = dp[1]
        B[2,i+1] = dp[0]
    elif self.rank == 3:
      for iNel,dp in enumerate( dphi ):
        i = 3 * iNel
        
        B[0,i  ] = dp[0]
        B[1,i+1] = dp[1]
        B[2,i+2] = dp[2]

        B[3,i+1] = dp[2]
        B[3,i+2] = dp[1]

        B[4,i  ] = dp[2]
        B[4,i+2] = dp[0]

        B[5,i  ] = dp[1]
        B[5,i+1] = dp[0]

    return B    

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def stress2matrix( self , stress ):

    T = zeros( shape=( self.rank*self.rank , self.rank*self.rank ) )

    if self.rank == 2:
      T[0,0] = stress[0]
      T[1,1] = stress[1]
      T[0,1] = stress[2]
      T[1,0] = stress[2]

      T[self.rank:,self.rank:] = T[:self.rank,:self.rank]
    elif self.rank == 3:
      T[0,0] = stress[0]
      T[1,1] = stress[1]
      T[2,2] = stress[2]
      T[1,2] = stress[3]
      T[0,2] = stress[4]
      T[0,1] = stress[5]
      T[2,1] = stress[3]
      T[2,0] = stress[4]
      T[1,0] = stress[5]

      T[self.rank:2*self.rank,self.rank:2*self.rank] = T[:self.rank,:self.rank]
      T[2*self.rank:,2*self.rank:] = T[:self.rank,:self.rank]
       
    return T

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBNLmatrix( self , dphi ):

    Bnl = zeros( shape=( self.rank*self.rank , self.rank*len(dphi) ) )

    if self.rank == 2:
      for i,dp in enumerate( dphi ):
        Bnl[0,2*i  ] = dp[0]
        Bnl[1,2*i  ] = dp[1]
        Bnl[2,2*i+1] = dp[0]
        Bnl[3,2*i+1] = dp[1]
    elif self.rank == 3:
      for i,dp in enumerate( dphi ):
        Bnl[0,3*i  ] = dp[0]
        Bnl[1,3*i  ] = dp[1]
        Bnl[2,3*i  ] = dp[2]

        Bnl[3,3*i+1] = dp[0]
        Bnl[4,3*i+1] = dp[1]
        Bnl[5,3*i+1] = dp[2]

        Bnl[6,3*i+2] = dp[0]
        Bnl[7,3*i+2] = dp[1]
        Bnl[8,3*i+2] = dp[2]

    return Bnl	    

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( self.rank , self.rank*len(h) ) )

    for i,a in enumerate( h ):
      for j in range(self.rank):
        N[j,self.rank*i+j] = a
    
    return N
