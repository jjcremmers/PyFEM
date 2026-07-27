# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseFailure import BaseFailure
from pyfem.materials.MatUtils    import vonMisesStress

class VonMises( BaseFailure ):

  def __init__ ( self, props ):
    """Initialize the Von Mises failure criterion.

    Parameters
    ----------
    props : object
      Material property container with the allowable stress ``smax``.
    """

    BaseFailure.__init__( self, props )
    
    self.smax = props.smax

  def check( self, stress , deformation ):
    """Compute the Von Mises failure index for the current stress state.

    Parameters
    ----------
    stress : array_like
      Stress components at the evaluation point.
    deformation : array_like
      Deformation measures for the current state. This argument is accepted
      for interface consistency and is not used in the calculation.

    Returns
    -------
    float
      Ratio of the Von Mises equivalent stress to the allowable stress.
    """

    FI = vonMisesStress( stress ) / self.smax
    
    return FI
