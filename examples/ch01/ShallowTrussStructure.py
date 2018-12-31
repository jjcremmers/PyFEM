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
#  Description: The Python file of the example presented in section        #
#               1.6 of the book, pages 25--29                              #
#                                                                          #
#  Usage:       python ShallowTrussStructure.py                            #
############################################################################

#############################
# Define problem parameters #
# see Figure 1.1            #
#############################

#Problem dimensions
b  = 10. #m
h  = 0.5 #m

#Material parameters
k   = 1000. #N/m
EA0 = 5e6   #N/m^2

#External load increment
DF = 50 #N

#Solver parameters
N       = 30
tol     = 1e-6
iterMax = 5

#Some useful functions
from math import sqrt

l    = lambda v : sqrt( b**2 +(h-v)**2 )
F    = lambda v : -EA0 * (h-v)/l(v) * (l(v)-l(0))/l(0) + k * v
dFdv = lambda v : (EA0/l(v)) * ((h-v)/l(v))**2  + k + (EA0/l(v)) * (l(v)-l(0))/l(0)

#############################
# Newton-Raphson iterations #
#############################

#Initialize
v    = 0.
Dv   = 0.
Fext = 0.

DFext = DF

output = [ [0.,0.] ]

#Load step iterator
for i in range(N):
  
  print('=================================')
  print(' Load step %i' % i)
  print('=================================')
  print('  NR iter : |Fext-F(v)|')

  #####################################  
  # Compute the new external force    #
  #####################################

  Fext = Fext + DFext
  
  #Initialize Newton-Raphson iteration parameters  
  error = 1.
  iiter = 0

  #############################################
  # Compute the derivative of F w.r.t v       #
  #############################################

  while error > tol:
    
    ########################################  
    # Solve for dv                         #
    ########################################

    dv = (1./dFdv(v+Dv))* ( Fext - F(v+Dv) )
    
    ############################################
    # Update Delta v                           #
    ############################################

    Dv += dv
    
    ###############################################
    # Convergence check                           #
    ###############################################

    error  = abs( Fext-F(v+Dv) ) 

    #Increment the Newton-Raphson iteration counter
    iiter += 1

    print('  Iter', iiter, ':', error)

    if iiter == iterMax:
      raise RuntimeError('Newton-Raphson iterations did not converge!')

  #Update the displacement
  v  += Dv
  Dv  = 0.

  #Store the output
  output.append( [ v, F(v) ] )

  print('=================================')
  

###############################
# Post-processing             #
###############################

from pylab import plot, show, xlabel, ylabel

plot( [x[0] for x in output], [x[1] for x in output], 'ro' )

#Exact solution
from numpy import arange

vrange = arange(0,1.2,0.01)
plot( vrange, [F(vval) for vval in vrange], 'b-' ) 
xlabel('v [m]')
ylabel('F [N]')

show()
