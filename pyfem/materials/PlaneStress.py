# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros, dot

class PlaneStress( BaseMaterial ):

  def __init__ ( self, props ):

    BaseMaterial.__init__( self, props )

    self.H = zeros( (3,3) )

    self.H[0,0] = self.E/(1.-self.nu*self.nu)
    self.H[0,1] = self.H[0,0]*self.nu
    self.H[1,0] = self.H[0,1]
    self.H[1,1] = self.H[0,0]
    self.H[2,2] = self.E/(2.0*(1.0+self.nu))

    self.outLabels = [ "S11" , "S22" , "S12" ]

  def getStress( self, deformation ):
    
    sigma = dot( self.H, deformation.strain )

    self.outData = sigma
    
    return sigma, self.H

  def getTangent( self ):
  
    return self.H

