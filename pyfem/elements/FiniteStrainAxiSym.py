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

from numpy import zeros, dot, outer, ones, eye, reshape
from math  import pi

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class FiniteStrainAxiSym( Element ):
  
  def __init__ ( self, elnodes , props ):
  
    self.method = "TL"
    
    Element.__init__( self, elnodes , props )

    if props.rank != 2:
      print("Error")
      
    self.dofTypes = [ 'u' , 'v' ]
    self.nstr     = 6
    
    self.kin = Kinematics( self.rank,self.nstr)
         
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
      print("Error")

#
#
#

  def getTLTangentStiffness ( self, elemdat ):
  
    sData = getElemShapeData( elemdat.coords )
   
    for iData in sData:

      r      = dot( elemdat.coords[:,0] , iData.h )
      weight = 2.0*pi*r*iData.weight
      
      kin = self.getKinematics( iData.dhdx , iData.h , elemdat , r ) 
      B   = self.getBmatrix   ( iData.dhdx , iData.h , kin.F , r )
      
      sigma,tang = self.mat.getStress( kin )
      
      t4 = self.tang6to4  ( tang  )
      s4 = self.stress6to4( sigma )
        
      stiff = dot ( B.transpose() , dot ( t4 , B ) )
      
      for i,dpi in enumerate(iData.dhdx):
        for j,dpj in enumerate(iData.dhdx):
          stiff[2*i,2*j] += s4[0]*dpi[0]*dpj[0] + s4[1]*dpi[1]*dpj[1]
          stiff[2*i,2*j] += s4[2]*iData.h[i]*iData.h[j]/(r*r)                 
          stiff[2*i,2*j] += s4[3]*(dpi[0]*dpj[1] + dpi[1]*dpj[0])
          
          stiff[2*i+1,2*j+1] += s4[0]*dpi[0]*dpj[0] + s4[1]*dpi[1]*dpj[1]      
          stiff[2*i+1,2*j+1] += s4[3]*(dpi[0]*dpj[1] + dpi[1]*dpj[0])
          
      elemdat.stiff += stiff * weight     
      elemdat.fint  += dot ( B.transpose() , s4 ) * weight
      
      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() )
      
#
#
#

  def getULTangentStiffness ( self, elemdat ):
  
    state0  = elemdat.state - elemdat.Dstate
    curCrds = elemdat.coords + reshape(state0,elemdat.coords.shape)
    
    sData0 = getElemShapeData( elemdat.coords )
    sDataC = getElemShapeData( curCrds )
   
    for iOrig,iCurr in zip(sData0,sDataC):

      r      = dot( curCrds[:,0] , iOrig.h )
      r0     = dot( elemdat.coords[:,0] , iOrig.h )
      
      weight = 2.0*pi*r*iCurr.weight
      
      kin = self.getKinematics( iOrig.dhdx , iOrig.h , elemdat , r0 ) 
      B   = self.getULBmatrix ( iCurr.dhdx , iCurr.h , r )
      
      sigma,tang = self.mat.getStress( kin )
      
      t4 = self.tang6to4  ( tang  )
      s4 = self.stress6to4( sigma )
        
      stiff = dot ( B.transpose() , dot ( t4 , B ) )
      
      for i,dpi in enumerate(iCurr.dhdx):
        for j,dpj in enumerate(iCurr.dhdx):
          stiff[2*i,2*j] += s4[0]*dpi[0]*dpj[0] + s4[1]*dpi[1]*dpj[1]
          stiff[2*i,2*j] += s4[2]*iCurr.h[i]*iCurr.h[j]/(r*r)                 
          stiff[2*i,2*j] += s4[3]*(dpi[0]*dpj[1] + dpi[1]*dpj[0])
          
          stiff[2*i+1,2*j+1] += s4[0]*dpi[0]*dpj[0] + s4[1]*dpi[1]*dpj[1]      
          stiff[2*i+1,2*j+1] += s4[3]*(dpi[0]*dpj[1] + dpi[1]*dpj[0])
          
      elemdat.stiff += stiff * weight     
      elemdat.fint  += dot ( B.transpose() , s4 ) * weight
      
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

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getTLInternalForce ( self, elemdat ):
   
    n = self.dofCount()
   
    sData = getElemShapeData( elemdat.coords )

    for iData in sData:
    
      r      = dot( elemdat.coords[:,0] , iData.h )
      weight = 2.0*pi*r*iData.weight
      
      kin = self.getKinematics( iData.dhdx , iData.h , elemdat , r ) 
      B   = self.getBmatrix   ( iData.dhdx , iData.h , kin.F , r )
                  
      sigma,tang = self.mat.getStress( kin )
       
      s4 = self.stress6to4( sigma )
      
      elemdat.fint    += dot ( B.transpose() , s4 ) * weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() ) 
      
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getULInternalForce ( self, elemdat ):
   
    n = self.dofCount()
    
    state0  = elemdat.state - elemdat.Dstate
    curCrds = elemdat.coords + reshape(state0,elemdat.coords.shape)
    
    sData0 = getElemShapeData( elemdat.coords )
    sDataC = getElemShapeData( curCrds )
   
    for iOrig,iCurr in zip(sData0,sDataC):
    
      r      = dot( curCrds[:,0] , iOrig.h )
      r0     = dot( elemdat.coords[:,0] , iOrig.h )
      
      weight = 2.0*pi*r*iCurr.weight
      
      kin = self.getKinematics( iOrig.dhdx , iOrig.h , elemdat , r0 ) 
      B   = self.getULBmatrix ( iCurr.dhdx , iCurr.h , r )
                  
      sigma,tang = self.mat.getStress( kin )
       
      s4 = self.stress6to4( sigma )
      
      elemdat.fint    += dot ( B.transpose() , s4 ) * weight

      self.appendNodalOutput( self.mat.outLabels() , self.mat.outData() ) 

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getMassMatrix ( self, elemdat ):
      
    sData = getElemShapeData( elemdat.coords )

    rho = elemdat.matprops.rho

    for iData in sData:
    
      r      = dot( elemdat.coords[:,0] , iData.h )
      weight = 2.0*pi*r*iData.weight
      
      N  = self.getNmatrix( iData.h )
      elemdat.mass += dot ( N.transpose() , N ) * rho * weight
     
    elemdat.lumped = sum(elemdat.mass)
            
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getKinematics( self , dphi , h , elemdat , r ):
  
    kin = Kinematics(3,6)

    kin.F  = eye(3)
    kin.F0 = eye(3)
    
    elstate  = elemdat.state
    elstate0 = elstate - elemdat.Dstate
    
    invr = 1.0/r
  
    for i,(dp,p) in enumerate(zip(dphi,h)):
      kin.F[0,0] += dp[0]*elstate[2*i]
      kin.F[0,1] += dp[1]*elstate[2*i]
      kin.F[1,0] += dp[0]*elstate[2*i+1]
      kin.F[1,1] += dp[1]*elstate[2*i+1]
      kin.F[2,2] += p*elstate[2*i] * invr
      
      kin.F0[0,0] += dp[0]*elstate0[2*i]
      kin.F0[0,1] += dp[1]*elstate0[2*i]
      kin.F0[1,0] += dp[0]*elstate0[2*i+1]
      kin.F0[1,1] += dp[1]*elstate0[2*i+1]
      kin.F0[2,2] += p*elstate0[2*i] * invr
      
    kin.E  = 0.5*(dot(kin.F.transpose(),kin.F)-eye(3))
    kin.E0 = 0.5*(dot(kin.F0.transpose(),kin.F0)-eye(3))

    dE = kin.E - kin.E0
    
    kin.strain[0] =     kin.E[0,0]
    kin.strain[1] =     kin.E[1,1]
    kin.strain[2] =     kin.E[2,2]
    kin.strain[3] = 2.0*kin.E[1,2]
    kin.strain[4] = 2.0*kin.E[0,2]
    kin.strain[5] = 2.0*kin.E[0,1]
    
    return kin

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getBmatrix( self , dphi , phi , F , r ):

    B = zeros( shape=(4, 2*len(dphi) ) )

    invr = 1.0/r
    
    for i,(dp,p) in enumerate(zip(dphi,phi)):
      B[0,2*i  ] = dp[0]*F[0,0]
      B[0,2*i+1] = dp[0]*F[1,0]
	
      B[1,2*i  ] = dp[1]*F[0,1]
      B[1,2*i+1] = dp[1]*F[1,1]
      
      B[2,2*i  ] = p * F[2,2] * invr
      B[2,2*i+1] = 0.0

      B[3,2*i  ] = dp[1]*F[0,0]+dp[0]*F[0,1]
      B[3,2*i+1] = dp[0]*F[1,1]+dp[1]*F[1,0]
 
    return B
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getULBmatrix( self , dphi , phi , r ):

    B = zeros( shape=(4, 2*len(dphi) ) )

    invr = 1.0/r
    
    for i,(dp,p) in enumerate(zip(dphi,phi)):
      B[0,2*i  ] = dp[0]
      B[1,2*i+1] = dp[1]
      
      B[2,2*i  ] = p * invr
      B[2,2*i+1] = 0.0

      B[3,2*i  ] = dp[1]
      B[3,2*i+1] = dp[0]
 
    return B
        

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getNmatrix( self , h ):

    N = zeros( shape=( 2 , 2*len(h) ) )

    for i,a in enumerate( h ):
      for j in range(2):
        N[j,2*i+j] = a
    
    return N
    
 #------------------------------------------------------------------------------
 #
 #------------------------------------------------------------------------------
    
  def stress6to4( self , sigma ):
  
    s4 = zeros(4)
    
    s4[0] = sigma[0]
    s4[1] = sigma[1]
    s4[2] = sigma[2]
    s4[3] = sigma[5]
    
    return s4
    
 #------------------------------------------------------------------------------
 #
 #------------------------------------------------------------------------------
 
  def tang6to4( self , tang ):
  
    t4 = zeros( shape=(4,4) )
    
    t4[:3,:3] = tang[:3,:3]
    t4[:3,3 ] = tang[:3,5]
    t4[3,:3]  = tang[5,:3]
    t4[3,3]   = tang[5,5]
    
    return t4
    
    
    
