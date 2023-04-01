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

from numpy import zeros,ones,dot,transpose
from numpy.linalg import inv
from math import sin,cos,pi,sqrt,tan,atan

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class TransverseIsotropic:
  
  def __init__( self , props ):

    if hasattr( props , "E1" ):
      self.E1 = props.E1
      self.E2 = props.E2
    elif hasattr( props , "E" ):
      if type(props.E) is list:
        if len(props.E) == 2:
          self.E1 = props.E[0]
          self.E2 = props.E[1]
        elif len(E) == 1:
          self.E1 = props.E[0]
          self.E2 = props.E[0]
        else:
          raise RuntimeError("Please add E as a float or a list; E = (E1,E2).")
      else:
        self.E1    = props.E
        self.E2    = props.E

    if hasattr( props , "nu12" ):
      self.nu12  = props.nu12
    else:
      self.nu12  = props.nu

    if hasattr( props , "G12" ):
      self.G12   = props.G12    
    else:
      self.G12   = self.E1/(2.0*(1.0+self.nu12))

    if hasattr( props , "G13" ):
      self.Q44   = props.G13    
    else:
      self.Q44   = self.G12

    if hasattr( props , "G23" ):
      self.Q55   = props.G23    
    else:
      self.Q55   = self.G12
  
    self.nu21  = self.E2/self.E1*self.nu12

    self.rho = props.rho

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getQ( self ):

    if not hasattr( self , 'Q' ):
      self.Q = zeros( shape=(3,3) )

      self.Q[0,0] = self.E1/(1.-self.nu12*self.nu21)
      self.Q[0,1] = self.nu12*self.E2/(1.0-self.nu12*self.nu21)
      self.Q[1,1] = self.E2/(1.-self.nu12*self.nu21)
      self.Q[1,0] = self.Q[0,1]
      self.Q[2,2] = self.G12
  
    return self.Q

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getU( self ):

    if not hasattr( self , 'U' ):
      self.getQ()

      self.U = zeros(5)

      self.U[0] = 0.125*(3.*self.Q[0,0]+3.*self.Q[1,1]+2.*self.Q[0,1]+4.*self.Q[2,2])
      self.U[1] = 0.5*(self.Q[0,0]-self.Q[1,1])
      self.U[2] = 0.125*(self.Q[0,0]+self.Q[1,1]-2.*self.Q[0,1]-4.*self.Q[2,2])
      self.U[3] = 0.125*(self.Q[0,0]+self.Q[1,1]+6.*self.Q[0,1]-4.*self.Q[2,2])
      self.U[4] = 0.5*(self.U[0]-self.U[3])

    return self.U

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getS( self ):

    self.S = zeros( shape=(3,3) )

    self.S[0,0] = 1./self.E1
    self.S[0,1] = -self.nu12/self.E1
    self.S[1,1] = 1./self.E2
    self.S[1,0] = self.S[0,1]
    self.S[2,2] = 1./self.G12

    return self.S

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
  
  def getV( self ):

    if not hasattr( self , 'V' ):
      self.getS()

      self.V = zeros(5)

      self.V[0] = 0.125*(3.*self.S[0,0]+3.*self.S[1,1]+2.*self.S[0,1]+self.S[2,2])
      self.V[1] = 0.5*(self.S[0,0]-self.S[1,1])
      self.V[2] = 0.125*(self.S[0,0]+self.S[1,1]-2.*self.S[0,1]-self.S[2,2])
      self.V[3] = 0.125*(self.S[0,0]+self.S[1,1]+6.*self.S[0,1]-self.S[2,2])
      self.V[4] = 2.*(self.V[0]-self.V[3])

    return self.V

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getQbar( self , theta ):

    if not hasattr( self , 'U' ):
      self.getU()

    Qbar = zeros( shape=(3,3) )

    rad = theta*pi/180.

    s2 = sin(2.*rad)
    s4 = sin(4.*rad)

    c2 = cos(2.*rad)
    c4 = cos(4.*rad)

    Qbar[0,0] = self.U[0]+self.U[1]*c2+self.U[2]*c4
    Qbar[0,1] = self.U[3]-self.U[2]*c4
    Qbar[1,0] = Qbar[0,1]
    Qbar[1,1] = self.U[0]-self.U[1]*c2+self.U[2]*c4
    Qbar[0,2] = 0.5*self.U[1]*s2+self.U[2]*s4
    Qbar[1,2] = 0.5*self.U[1]*s2-self.U[2]*s4
    Qbar[2,0] = Qbar[0,2]
    Qbar[2,1] = Qbar[1,2]
    Qbar[2,2] = self.U[4]-self.U[2]*c4

    return Qbar

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getQshearbar( self , theta ):

    Qshear = zeros( shape=(2,2) )
   
    rad = theta*pi/180.

    Qshear[0,0] = self.Q44*cos(rad)*cos(rad)+self.Q55*sin(rad)*sin(rad)
    Qshear[1,1] = self.Q55*cos(rad)*cos(rad)+self.Q44*sin(rad)*sin(rad)
    Qshear[0,1] = (self.Q55-self.Q44)*cos(rad)*sin(rad)
    Qshear[1,0] = Qshear[0,1]

    return Qshear

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getSbar( self , theta ):

    if not hasattr( self , 'V' ):
      self.getV()

    Sbar = zeros( shape=(3,3) )

    rad = theta*pi/180.

    s2 = sin(2.*rad)
    s4 = sin(4.*rad)

    c2 = cos(2.*rad)
    c4 = cos(4.*rad)

    Sbar[0,0] = self.V[0]+self.V[1]*c2+self.V[2]*c4
    Sbar[0,1] = self.V[3]-self.V[2]*c4
    Sbar[1,0] = Sbar[0,1]
    Sbar[1,1] = self.V[0]-self.V[1]*c2+self.V[2]*c4
    Sbar[0,2] = self.V[1]*s2+2.*self.V[2]*s4
    Sbar[1,2] = self.V[1]*s2-2.*self.V[2]*s4
    Sbar[2,0] = Sbar[0,2]
    Sbar[2,1] = Sbar[1,2]
    Sbar[2,2] = self.V[4]-4.*self.V[2]*c4

    return Sbar

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class Layer:

  def __init__( self , props ):
   
    if type(props.material) == str:
      self.mat   = props.material
    else:
      self.mat   = "mat"

    if hasattr( props , "theta" ):
      self.theta = props.theta    
    else:
      self.theta = 0.0
    self.thick = props.thickness

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class Laminate:

  def __init__( self , props ): 

    self.materials = {}
    self.layers    = []

    if hasattr( props , "layers" ):

      matNames = props.materials
      layNames = props.layers

      for name in matNames:
        self.materials[name] = TransverseIsotropic( getattr( props, name ) )

      for layer in layNames:
        self.layers.append( Layer( getattr( props, layer ) ) )
        
    else:
      self.materials["mat"] = TransverseIsotropic( props.material )
      self.layers.append( Layer( props ) )
    
    self.h     = zeros( len(self.layers)+1 )
    self.thick = 0.
    
    for i,layer in enumerate(self.layers):
      self.h[i+1] = self.thick+layer.thick
      self.thick += layer.thick

    self.h += -0.5*self.thick*ones( len(self.h) )

    self.shearCorr = 5.0/6.0

    if hasattr( props , "shearCorrection" ):
      self.shearCorr = props.shearCorrection

#------------------------------------------------------------------------------
#  
#------------------------------------------------------------------------------

  def layerCount( self ):
   
    return len(self.layers)

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getA( self ):

    self.A = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.mat
      theta = layer.theta
      
      self.A 	+= self.materials[name].getQbar( theta ) * (self.h[i+1]-self.h[i])

    return self.A

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getB( self ):

    self.B = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.mat
      theta = layer.theta

      self.B += 0.5*self.materials[name].getQbar( theta ) * (self.h[i+1]**2-self.h[i]**2)

    return self.B

  def getD( self ):

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    self.D = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.mat
      theta = layer.theta

      self.D += 1.0/3.0*self.materials[name].getQbar( theta ) * (self.h[i+1]**3-self.h[i]**3)

    return self.D

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getQbar( self , i ):

    name  = self.layers[i].mat
    theta = self.layers[i].theta

    return self.materials[name].getQbar( theta )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getQ( self , i ):

    name  = self.layers[i].mat

    return self.materials[name].getQ()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getAshear( self ):
  
    self.Ashear = zeros( shape=(2,2) )
    
    for i,layer in enumerate(self.layers):
      name  = layer.mat
      theta = layer.theta

      self.Ashear += self.shearCorr*self.materials[name].getQshearbar( theta )*(self.h[i+1]-self.h[i])
    
    return self.Ashear

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getMassInertia( self ):

    massInert = zeros( 3 )

    for i,layer in enumerate(self.layers):
      name  = layer.mat

      massInert[0] += self.materials[name].rho*(self.h[i+1]-self.h[i])
      massInert[1] += 0.5*self.materials[name].rho*(self.h[i+1]**2-self.h[i]**2)
      massInert[2] += 1./3.*self.materials[name].rho*(self.h[i+1]**3-self.h[i]**3)
    
    return massInert

#==============================================================================
#  Utility functions
#==============================================================================

#------------------------------------------------------------------------------
#  stressTransformation 
#    Transforms stress from 12 coordinate system to xy coordinate system
#------------------------------------------------------------------------------

def stressTransformation( sigma , theta ):

  signew = zeros( 3 )

  rad = theta*pi/180.

  c = cos(rad)
  s = sin(rad)

  signew[0] = sigma[0]*c*c + sigma[1]*s*s + 2.*sigma[2]*c*s
  signew[1] = sigma[0]*s*s + sigma[1]*c*c - 2.*sigma[2]*c*s
  signew[2] = ( sigma[1] - sigma[0] )*c*s + sigma[2]*( c*c - s*s )
   
  return signew


