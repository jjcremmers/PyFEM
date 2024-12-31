################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
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

