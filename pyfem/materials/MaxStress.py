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

    if len(stress) = 3:  
      FI = max( fabs( stress[2] / self.S ) , FI )
    else:
      FI = max( fabs( stress[5] / self.S ) , FI )

    return FI
