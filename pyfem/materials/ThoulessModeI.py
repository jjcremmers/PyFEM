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

    self.outData = trac

    return stress,tang
