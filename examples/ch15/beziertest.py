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
############################################################################
#  Description: The Python file of the bezier interpolation example as     #
#               presented in section 15.3, page 487-490.                   #
#                                                                          #
#  Usage:       python beziertest.py                                       #
############################################################################

from numpy import zeros,dot

C = zeros(shape=[4,4,4])
   
C[0,0,0] = 1.0
C[0,1,1] = 1.0
C[0,1,2] = 0.5
C[0,1,3] = 0.25
C[0,2,2] = 0.5
C[0,2,3] = 0.5
C[0,3,3] = 0.25

C[1,0,0] = 0.25
C[1,1,0] = 0.5
C[1,1,1] = 0.5
C[1,2,0] = 0.25
C[1,2,1] = 0.5
C[1,2,2] = 1.0
C[1,2,3] = 0.5
C[1,3,3] = 0.5

C[2,0,0] = 0.5
C[2,1,0] = 0.5
C[2,1,1] = 1.0
C[2,2,2] = 1.0
C[2,3,3] = 1.0

C[3,0,0] = 1.0
C[3,1,1] = 1.0
C[3,2,2] = 1.0
C[3,3,3] = 1.0

coords = zeros( shape=(10,2) )

coords[0,:] = [0.0 , 0.0 ]
coords[1,:] = [1.0 , 1.0 ]
coords[2,:] = [1.0 , 3.0 ]
coords[3,:] = [2.0 , 3.0 ]
coords[4,:] = [2.5 , 1.5 ]
coords[5,:] = [1.5 , 0.5 ]
coords[6,:] = [3.0 , 0.0 ]
coords[7,:] = [3.2 , 2.0 ]
coords[8,:] = [3.8 , 2.5 ]
coords[9,:] = [4.0 , 0.0 ]

elems = zeros( shape=(4,4) , dtype=int )

elems[0,:] = [ 0 , 1 , 2 , 3 ]
elems[1,:] = [ 1 , 2 , 3 , 4 ]
elems[2,:] = [ 3 , 4 , 5 , 6 ]
elems[3,:] = [ 6 , 7 , 8 , 9 ]

output = []
length = 0.

from pyfem.util.BezierShapeFunctions import getElemBezierData

for elemNodes,Celem in zip(elems,C):

  sdata = getElemBezierData ( coords[elemNodes,:] , Celem , \
                              order = 100 , elemType='Line4')
 
  for idata in sdata:
    
    x = dot(idata.h,coords[elemNodes,:])

    output.append(x)

    length += idata.weight

print("The length of the curve is ",length)

from pylab import plot, show, xlabel, ylabel

plot( [x[0] for x in output], [x[1] for x in output], '-' )
plot( [x[0] for x in coords], [x[1] for x in coords], 'ro-' )
 
show()
    
