# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.materials.BaseMaterial import BaseMaterial
from pyfem.materials.MaterialManager import MaterialManager

from numpy import zeros, dot

class MultiMaterial( BaseMaterial ):

  def __init__ ( self, props ):

    BaseMaterial.__init__( self, props )
 
    self.matmodels = []
    
    for material in self.materials:    
      matProps            = getattr( props , material )
      matProps.rank       = props.rank
      matProps.solverStat = self.solverStat
  
      matmodel = MaterialManager( getattr( props , material ) )
      
      self.matmodels.append( matmodel )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getStress( self, deformation ):

    stress, tang = self.matmodels[deformation.iMat].getStress( deformation )
    
    self.outData = self.matmodels[deformation.iMat].outData()
    
    return stress , tang
