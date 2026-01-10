# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseFailure import BaseFailure
from pyfem.materials.MatUtils    import vonMisesStress

class VonMises( BaseFailure ):

  def __init__ ( self, props ):

    BaseFailure.__init__( self, props )
    
    self.smax = props.smax

  def check( self, stress , deformation ):

    FI = vonMisesStress( stress ) / self.smax
    
    print(FI)
    return FI
