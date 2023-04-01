############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by            #
#  Joris J.C. Remmers, Clemens V. Verhoosel and Rene de Borst and since    #
#  then augmented and  maintained by Joris J.C. Remmers.                   #
#  All rights reserved.                                                    #
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
from math import pi

import sys

class ThermoSmallStrainAxiSym( Element ):
  
  def __init__ ( self, elnodes , props ):

    Element.__init__( self, elnodes , props )

    if props.rank != 2: 
      raise RuntimeError("This is an axisymmetric element. Please use an input mesh with rank 2.") 
      
    self.dofTypes = [ 'u' , 'v' , 'temp' ]
  
    self.kin = Kinematics(3,6)
    
    self.D     = self.material.heatConductivity*eye(2)
    self.capac = self.material.heatCapacity
    
    self.alpha = self.material.alpha*ones(4)
    self.alpha[3] = 0.;
    
    self.labels = [ "qr" , "qz" ]
    
    self.transient = True
    self.theta = 1.0
      
  def __type__ ( self ):
    return name
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):
       
    sData = getElemShapeData( elemdat.coords )
    
    dDofs,tDofs = self.splitDofIDs( len(elemdat.coords) )
    
    temp0 = elemdat.state [tDofs] - elemdat.Dstate[tDofs]
    
    if self.transient:
      ctt      = zeros(shape=(len(tDofs),len(tDofs)))
      invdtime = 1.0/self.solverStat.dtime
                       
    for iInt,iData in enumerate(sData):
      
      r      = dot( elemdat.coords[:,0] , iData.h )
      weight = 2.0*pi*r*iData.weight
      
      B = self.getBmatrix( iData.dhdx , iData.h , r )

      strain  = dot ( B , elemdat.state [dDofs] )
      dstrain = dot ( B , elemdat.Dstate[dDofs] )
      
      self.kin.strain  = self.axisymTo3D( strain  )
      self.kin.dstrain = self.axisymTo3D( dstrain ) 
      
      temp     = sum( iData.h * elemdat.state [tDofs] )
      dtemp    = sum( iData.h * elemdat.Dstate[tDofs] )
      gradTemp = dot( iData.dhdx.transpose() , elemdat.state [tDofs] )
            
      self.kin.strain[:3]  += -self.alpha[:3] * temp
      self.kin.dstrain[:3] += -self.alpha[:3] * dtemp 
            
      sigma,tang = self.mat.getStress( self.kin )
      
      s4 = self.stress6to4( sigma )
      t4 = self.tang6to4  ( tang )
      
      elemdat.stiff[ix_(dDofs,dDofs)] += \
        dot ( B.transpose() , dot ( t4 , B ) ) * weight
              
      elemdat.stiff[ix_(dDofs,tDofs)] += \
        dot ( B.transpose() , outer ( -self.alpha , iData.h ) ) * weight
      
      elemdat.stiff[ix_(tDofs,tDofs)] += \
        dot ( iData.dhdx , dot( self.D , iData.dhdx.transpose() ) ) * weight
  
      elemdat.fint[dDofs] += dot ( B.transpose() , s4 ) * weight
      
      if self.transient:
        ctt += self.capac * outer( iData.h , iData.h ) * weight
              
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
      self.appendNodalOutput( self.labels , dot(self.D,gradTemp) ) 
    
    if self.transient:  
      ktt0 = invdtime * ctt - elemdat.stiff[ix_(tDofs,tDofs)] * \
        ( 1.0-self.theta )
      
      elemdat.stiff *= self.theta
      
      elemdat.stiff[ix_(tDofs,tDofs)] += invdtime * ctt 
        
    elemdat.fint[tDofs] += \
      dot ( elemdat.stiff[ix_(tDofs,tDofs)] , elemdat.state[tDofs] )
      
    if self.transient:
      elemdat.fint[tDofs] += -dot ( ktt0 , temp0 )
     
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):
     
    sData = getElemShapeData( elemdat.coords )
    
    dDofs,tDofs = self.splitDofIDs( len(elemdat.coords) )
    
    temp0 = elemdat.state [tDofs] - elemdat.Dstate[tDofs]
    
    stiff = zeros(shape=(len(tDofs),len(tDofs)))
    
    if self.transient:
      ctt = zeros(shape=(len(tDofs),len(tDofs)))
      invdtime = 1.0/self.solverStat.dtime
                 
    for iInt,iData in enumerate(sData):
      
      B = self.getBmatrix( iData.dhdx )

      self.kin.strain  = dot ( B , elemdat.state [dDofs] )
      self.kin.dstrain = dot ( B , elemdat.Dstate[dDofs] )
      
      temp     = sum( iData.h * elemdat.state [tDofs] )
      dtemp    = sum( iData.h * elemdat.Dstate[tDofs] )
      gradTemp = dot( iData.dhdx.transpose() , elemdat.state [tDofs] )
            
      self.kin.strain[:3]  += -self.alpha[:3] * temp
      self.kin.dstrain[:3] += -self.alpha[:3] * dtemp 
            
      sigma,tang = self.mat.getStress( self.kin )
            
      stiff[ix_(tDofs,tDofs)] += \
        dot ( iData.dhdx , dot( self.D , iData.dhdx.transpose() ) ) * iData.weight
  
      elemdat.fint[dDofs] += dot ( B.transpose() , sigma ) * iData.weight
      
      if self.transient:
        ctt += self.capac * outer( iData.h , iData.h ) * iData.weight
              
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
      self.appendNodalOutput( self.labels , dot(self.D,gradTemp) ) 
    
    if self.transient:  
      stiff *= self.theta
      
      stiff += invdtime * ctt 
        
    elemdat.fint[tDofs] += dot ( stiff , elemdat.state[tDofs] )
      
    if self.transient:
      elemdat.fint[tDofs] += -dot ( ktt0 , temp0 )
       
#-------------------------------------------------------------------------------
#  getBmatrix
#-------------------------------------------------------------------------------

  def getBmatrix( self , dphi , phi , r ):

    b = zeros( shape=( 4 , 2*len(phi) ) )

    invr = 1.0/r
    
    for i,(dp,p) in enumerate(zip(dphi,phi)):
      b[0,i*2+0] = dp[0]
      b[1,i*2+1] = dp[1]
      b[2,i*2+0] = p*invr
      b[3,i*2+0] = dp[1]
      b[3,i*2+1] = dp[0]
 
    return b
    
#-----------------------------
#
#------------------------------

  def axisymTo3D( self , strain ):
  
    s6 = zeros(6)
    
    s6[0] = strain[0]
    s6[1] = strain[1]
    s6[2] = strain[2]
    s6[5] = strain[3]
    
    return s6

#----------------------------------------
#
#---------------------------------------

  def stress6to4( self , sigma ):
  
    s4 = zeros(4)
    
    s4[0] = sigma[0]
    s4[1] = sigma[1]
    s4[2] = sigma[2]
    s4[3] = sigma[5]
    
    return s4

#----------------------------------------
#
#---------------------------------------
    
  def tang6to4( self , tang ):
  
    t4 = zeros( shape=(4,4) )
    
    t4[:3,:3] = tang[:3,:3]
    t4[:3,3 ] = tang[:3,5]
    t4[3,:3]  = tang[5,:3]
    t4[3,3]   = tang[5,5]
    
    return t4

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def splitDofIDs( self , n ):
  
    '''Routine to split the dof IDs in two groups, one for the displacement 
       degrees of freedom, the second for the phase field degres of freedom. 
       n is the number of degrees of freedom in this model'''
    
    if n == 3:
      return [0,1,3,4,6,7],[2,5,8]
    elif n == 4:
      return [0,1,3,4,6,7,9,10],[2,5,8,11]
    else:
      print("Error")
