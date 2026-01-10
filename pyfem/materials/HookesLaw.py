# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseMaterial import BaseMaterial
import numpy as np

class HookesLaw( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    #Create the tangent matrix
    self.H = np.zeros( (1,1) )

    self.H(0,0) = self.E

  def getStress( self, deformation ):

    sigma = self.E * deformation.eps

    return sigma, self.H

  def getTangent( self , deformation ):
  
    return self.H

