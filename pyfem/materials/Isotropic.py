############################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:  #
#                                                                          #
#    'Non-Linear Finite Element Analysis of Solids and Structures'         #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel        #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                        #
#                                                                          #
#  The code is written by J.J.C. Remmers, C.V. Verhoosel and R. de Borst.  #
#                                                                          #
#  The latest stable version can be downloaded from the web-site:          #
#     http://www.wiley.com/go/deborst                                      #
#                                                                          #
#  A github repository, with the most up to date version of the code,      #
#  can be found here:                                                      #
#     https://github.com/jjcremmers/PyFEM                                  #
#                                                                          #
#  The code is open source and intended for educational and scientific     #
#  purposes only. If you use PyFEM in your research, the developers would  #
#  be grateful if you could cite the book.                                 #  
#                                                                          #
#  Disclaimer:                                                             #
#  The authors reserve all rights but do not guarantee that the code is    #
#  free from errors. Furthermore, the authors shall not be liable in any   #
#  event caused by the use of the program.                                 #
############################################################################

from pyfem.materials.BaseMaterial import BaseMaterial
from numpy import zeros, dot

class Isotropic( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    #Create the hookean matrix
    self.H = zeros( (6,6) )

    fac = 1.0 / (2.0 * self.nu * self.nu + self.nu - 1.0 );
  
    self.H[0,0] = fac * self.E * ( self.nu - 1.0 );
    self.H[0,1] = -1.0 * fac * self.E * self.nu;
    self.H[0,2] = self.H[0,1];
    self.H[1,0] = self.H[0,1];                                  
    self.H[1,1] = self.H[0,0];                                  
    self.H[1,2] = self.H[0,1];                                  
    self.H[2,0] = self.H[0,1];                                  
    self.H[2,1] = self.H[0,1];                                  
    self.H[2,2] = self.H[0,0];                                  
    self.H[3,3] = self.E / ( 2.0 + 2.0 * self.nu );          
    self.H[4,4] = self.H[3,3];                                  
    self.H[5,5] = self.H[3,3];

  def getStress( self, deformation ):

    sigma = dot( self.H, deformation.strain )

    return sigma, self.H

  def getTangent( self ):
  
    return self.H

