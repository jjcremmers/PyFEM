# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros
from math  import exp 

class ThoulessModeI( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    self.d3    = 2.*self.Gc / ( ( -self.d1d3+self.d2d3+1.0)*self.Tult )
    self.d1    = self.d1d3 * self.d3
    self.d2    = self.d2d3 * self.d3
    self.dummy = self.Tult / self.d1

    #Set the labels for the output data in this material model
    self.outLabels = [ "Tn" , "Ts" ]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getStress( self, deformation ):
 
    stress = zeros( 2 )
    tang   = zeros( (2,2) )

    if deformation.strain[0] < self.d1:
      stress[0] = self.dummy * deformation.strain[0]
      tang[0,0] = self.dummy
    elif deformation.strain[0] >= self.d1 and deformation.strain[0] < self.d2:
      stress[0] = self.Tult
      tang[0,0] = 0.0
    elif deformation.strain[0] >= self.d2 and deformation.strain[0] < self.d3:
      stress[0] = self.Tult* ( 1.0 - ( deformation.strain[0] - self.d2 ) / ( self.d3-self.d2 ) )
      tang[0,0] = self.Tult* ( - 1.0 ) / ( self.d3 - self.d2 )
    else:	
      stress[0] = 0.0
      tang[0,0] = 0.0

    self.outData = stress

    return stress,tang
