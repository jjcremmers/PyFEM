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
from numpy import zeros
from math  import exp 

class XuNeedleman( BaseMaterial ):

  def __init__ ( self, props ):

    #Call the BaseMaterial constructor
    BaseMaterial.__init__( self, props )

    self.r = 0.
    self.q = 1.

    self.vnmax = self.Gc/(2.71828183*self.Tult);
    self.vtmax = self.q*self.Gc/(1.16580058*self.Tult);

    #Set the labels for the output data in this material model
    self.outLabels = [ "Tn" , "Ts" ]

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def getStress( self, deformation ):
 
    stress = zeros( 2 )
    tang   = zeros( (2,2) )

    t1 = 1.0/self.vnmax
    t3 = deformation.strain[0]*t1
    t4 = exp(-t3)
    t6 = 1.0-self.q
    t9 = 1.0/(self.r-1.0)
    t12 = (self.r-self.q)*t9
    t14 = self.q+t12*t3
    t15 = deformation.strain[1]*deformation.strain[1]
    t16 = self.vtmax*self.vtmax
    t17 = 1.0/t16
    t18 = 0.0
    t19 = exp(-t15*t17)
    t24 = self.Gc*t4

    stress[0] = -t4*((1.0-self.r+t3)*t6*t9-t14*t19)*self.Gc*t1+t24*(t1*t6*t9-t12*t1*t19)
    stress[1] = 2.0*t24*t14*deformation.strain[1]*t17*t19

    t1  = self.vnmax*self.vnmax
    t4  = 1/self.vnmax
    t5  = deformation.strain[0]*t4
    t6  = exp(-t5)
    t8  = 1.0-self.q
    t11 = 1/(self.r-1.0)
    t14 = (self.r-self.q)*t11
    t16 = self.q+t14*t5
    t17 = deformation.strain[1]*deformation.strain[1]
    t18 = self.vtmax*self.vtmax
    t19 = 1/t18
    t21 = exp(-t17*t19)
    t26 = self.Gc*t4
    t38 = t19*t21
    t41 = self.Gc*t6
    t46 = -t26*t6*t16*deformation.strain[1]*t38+t41*t14*t4*deformation.strain[1]*t38
    t52 = t18*t18
    
    tang[0,0] = self.Gc/t1*t6*((1.0-self.r+t5)*t8*t11-t16*t21)-2.0*t26*t6*(t4*t8*t11-t14*t4*t21)
    tang[0,1] = 2.0*t46
    tang[1,0] = 2.0*t46
    tang[1,1] = 2.0*t41*t16*t19*t21-4.0*t41*t16*t17/t52*t21

    self.outData = stress

    return stress,tang
