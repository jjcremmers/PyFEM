# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros
from math  import exp 

class PowerLawModeI( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    self.deltan  = self.Gc / ( exp(1.0) * self.Tult )
    self.deltan2 = self.deltan *self.deltan
    self.deltan3 = self.deltan2*self.deltan

    self.outLabels = [ "Tn" , "Ts" ]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getStress( self, deformation ):
 
   stress = zeros( 2 )
   tang   = zeros( (2,2) )

   jump = deformation.strain[0]

   stress[0] = self.Gc/self.deltan2*exp(-jump/self.deltan)*jump
   tang[0,0] = self.Gc/self.deltan2*exp(-jump/self.deltan)*(1.0-jump/self.deltan)

   self.outData = stress

   return stress,tang
