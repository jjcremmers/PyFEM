# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from .Element import Element

from numpy import zeros, dot, eye, array
from scipy.linalg import norm
from math import atan2, sin, cos, tan

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class BeamLin ( Element ):

  #dofs per element
  dofTypes = [ 'u' , 'v' , 'w' , 'rz' , 'ry' , 'rz' ]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def __init__ ( self, elnodes , props ):
    Element.__init__( self, elnodes , props )

    self.EA  = self.E * self.A
    self.EI3 = self.E * self.I3
    self.EI2 = self.E * self.I2
    
    self.GA  = self.G * self.A
    self.GJ  = self.G * self.J
    
    self.family = "BEAM"
    
    self.bodyForce = False
    
    if hasattr(props,"bodyForce"):
      if props.bodyForce:
        self.bodyForce = True
        
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getTangentStiffness ( self, elemdat ):

    length , T = self.getT( elemdat )
        
    l2 = length*length   
    l3 = length*l2
    
    stiff = zeros( shape=(12,12))
    
    stiff[ 0, 0] = self.EA/length
    stiff[ 0, 6] = -stiff[0,0]
    
    stiff[ 1, 1] = 12*self.EJ3/l3
    stiff[ 1, 5] = 6*self.EJ3/l2   
    stiff[ 1, 7] = -stiff[1,1]
    stiff[ 1,11] =  stiff[1,5]
    
    stiff[ 2, 2] = 12*self.EJ2/l3
    stiff[ 2, 4] = -6*self.EJ2/l2
    stiff[ 2, 8] = -stiff[2,2]
    stiff[ 2,10] =  stiff[2,4]    
    
    stiff[ 3, 3] = self.GJ/length
    stiff[ 3, 9] = -stiff[3,3]
    
    stiff[ 4, 2] = -6*self.EJ2/l2
    stiff[ 4, 4] =  4*self.EJ2/length    
    stiff[ 4, 8] = -stiff[4,2]
    stiff[ 4,10] =  0.5*stiff[4,4]    
    
    stiff[ 5, 1] =  6*self.EJ3/l2
    stiff[ 5, 5] =  4*self.EJ3/length
    stiff[ 5, 7] = -stiff[5,1]
    stiff[ 5,11] =  0.5*stiff[5,5]    
    
    stiff[ 6, 0] = -stiff[0,0]
    stiff[ 6, 6] =  stiff[0,0]  
    
    stiff[ 7, 1] = -stiff[1,1]
    stiff[ 7, 5] = -stiff[1,5]
    stiff[ 7, 7] =  stiff[1,1]
    stiff[ 7,11] = -stiff[1,5]
    
    stiff[ 8, 2] = -stiff[2,2]
    stiff[ 8, 4] = -6*self.EJ2/l3    
    stiff[ 8, 8] = -stiff[8,2]
    stiff[ 8,10] =  stiff[8,4]    
    
    stiff[ 9, 3] = self.GJ/length
    stiff[ 9, 9] = -stiff[9,3]
    
    stiff[10, 2] = 12*self.EJ2/l3
    stiff[10, 4] = -6*self.EJ2/l3    
    stiff[10, 8] = -stiff[4,2]
    stiff[10,10] =  stiff[4,4]    
    
    stiff[11, 2] = 12*self.EJ2/l3
    stiff[11, 4] = -6*self.EJ2/l3    
    stiff[11, 8] = -stiff[5,2]
    stiff[11,10] =  stiff[5,4]          
    
    elemdat.stiff = self.loc2glob( stiff , T )       
    
    
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getInternalForce ( self, elemdat ):

    length , T = self.getT( elemdat )
        
    state = self.glob2loc( elemdat.state , T )
    
    f    = self.getF   ( state , length )
    fvar = self.getFvar( f , length )
    
    ae , kt1 = self.getTransformation( state , fvar , length ) 
    
    fint = dot( ae.transpose() , fvar )
                              
    elemdat.fint  = self.loc2glob( fint  , T )
    
#----------------------------------------------------------------------
    
  def getMassMatrix ( self, elemdat ):
      
    mass = zeros( shape=(6,6))

    length , T = self.getT( elemdat )
    
    mass[0,0] = 140.0
    mass[0,3] =  70.0    
    
    mass[1,1] = 156.0
    mass[1,2] =  22.0*length
    mass[1,4] =  54.0
    mass[1,5] = -13.0*length
    
    mass[2,2] =   4.0*length*length
    mass[2,4] =  13.0*length
    mass[2,5] =  -3.0*length*length
    
    mass[3,3] = 140.0
    
    mass[4,4] = 156.0
    mass[4,5] = -22.0*length
    
    mass[5,5] =   4.0*length*length
    
    mass[2,1] = mass[1,2]
    mass[3,0] = mass[0,3]
    
    mass[4,1] = mass[1,4]
    mass[4,2] = mass[2,4]
    
    mass[5,1] = mass[1,5]
    mass[5,2] = mass[2,5]
    mass[5,4] = mass[4,5]        
    
    mass *= self.rho*self.A*length/420.0
    
    elemdat.mass = self.loc2glob( mass , T )                  
    elemdat.lumped = sum(elemdat.mass)   
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def getExternalForce( self, elemdat ):
              
    if self.bodyForce:   

      length , T = self.getT( elemdat )
         
      g = array([0.0,-9.81])
      
      b = 0.5 * g * self.rho * self.A * length
      
      elemdat.fint[:2]  += b
      elemdat.fint[3:5] += b
                     
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
    
  def getT ( self, elemdat ):

    T = zeros(shape=(2,2))
    
    lvec = elemdat.coords[1,:] - elemdat.coords[0,:]
  
    length = norm(lvec)

    lvec   = lvec*1.0/length
              
    T[0,0] =  lvec[0]
    T[0,1] =  lvec[1]
    T[1,0] = -lvec[1]
    T[1,1] =  lvec[0]
    
    return length, T
  


#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def glob2loc( self , a , t ):

    b = zeros(6)

    b[2] = a[2];
    b[5] = a[5];

    b[0] = t[0,0]*a[0] + t[0,1] * a[1];
    b[1] = a[0]*t[1,0] + t[1,1] * a[1];

    b[3] = a[3]*t[0,0] + t[0,1] * a[4];
    b[4] = a[3]*t[1,0] + t[1,1] * a[4]; 
    
    return b

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def loc2glob( self , a , T ):

    tt = eye(6)

    tt[0:2,0:2] = T
    tt[3:5,3:5] = T
    
    if a.ndim == 1:
      return dot( tt.transpose() , a )
    else:
      return dot( tt.transpose() , dot( a , tt ) )
