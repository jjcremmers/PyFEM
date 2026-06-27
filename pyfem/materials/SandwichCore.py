# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros, dot

class SandwichCore( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    if not hasattr( props , "factor" ):
      self.factor = 0.001

    if not hasattr( props , "G13" ):
      self.G13 = self.G

    if not hasattr( props , "G23" ):
      self.G23 = self.G

    self.H = zeros( (6,6) )

    self.H[0,0] = self.factor * self.E3
    self.H[1,1] = self.factor * self.E3
    self.H[2,2] = self.E3
    self.H[3,3] = self.factor * (self.G13 + self. G23)*0.5                               
    self.H[4,4] = self.G23                                  
    self.H[5,5] = self.G13    

    #Set the labels for the output data in this material model
    self.outLabels = [ "S11" , "S22" , "S33" , "S23" , "S13" , "S12" ]                             

  def getStress( self, deformation ):

    sigma = dot( self.H, deformation.eps )

    self.outData = sigma

    return sigma, self.H

  def getTangent( self ):
  
    return self.H

