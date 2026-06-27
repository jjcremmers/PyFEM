# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseFailure import BaseFailure

class MaxStress( BaseFailure ):

  def __init__ ( self, props ):

    BaseFailure.__init__( self, props )
    
    self.Xt = props.Xt
    self.Xc = props.Xc
    self.Yt = props.Yt
    self.Yc = props.Yc
    self.S  = props.S
    
  def check( self, stress , deformation ):

    FI = 0.
    
    if stress[0] > 0.:
      FI = stress[0]/self.Xt
    else:
      FI = stress[0]/-self.Xc
   
    if stress[1] > 0.:
      FI = max( stress[1]/self.Yt , FI )
    else:
      FI = max( stress[1]/-self.Yc , FI )

    if len(stress) == 3:  
      FI = max( abs( stress[2] / self.S ) , FI )
    else:
      FI = max( abs( stress[5] / self.S ) , FI )

    return FI
