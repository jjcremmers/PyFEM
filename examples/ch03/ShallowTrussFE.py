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
#  Description: The Python file of the Finite Element implementation       #
#               of the shallow truss problem as presented in section 3.2,  #  
#               page 76--84.                                               #
#                                                                          #
#  Use:         python3 ShallowTrussFE.py                                  #
############################################################################

import sys,os
sys.path.insert(0, r"../.." )

#############################
# Define problem parameters #
# see Figure 1.1            #
#############################

E    = 5e6
Area = 1.0
k    = 1000.

b    = 10.0
h    = 0.5

#External load increment
DF   = 50

#Solver parameters
N       = 30
tol     = 1e-6
iterMax = 5

##############################
# Store data in properties   #
##############################

from pyfem.util.dataStructures import Properties

props = Properties()
props.TrussElem  = Properties( { 'type' : 'Truss'  , 'E' : 5e6 , 'Area' : 1.0 } )
props.SpringElem = Properties( { 'type' : 'Spring' , 'k' : 2000. } )

#############################
# Defining finite element   #
# data structures           #
#############################

#NodeSet
from pyfem.fem.NodeSet import NodeSet

nodes = NodeSet()

nodes.add( 1, [  0.,0.  ] ) #Attachement point of spring (can be positioned anywhere along the y-axis)
nodes.add( 2, [-10.,0.  ] ) #Left support
nodes.add( 3, [ 10.,0.  ] ) #Right support
nodes.add( 4, [ 0. ,0.5 ] ) #Loading point

#ElementSet
from pyfem.fem.ElementSet import ElementSet

elements = ElementSet( nodes , props )

elements.add( 1, 'TrussElem'  , [2,4]  )
elements.add( 2, 'TrussElem'  , [3,4]  )
elements.add( 3, 'SpringElem' , [1,4]  )

#DofSpace
from pyfem.fem.DofSpace import DofSpace

dofs = DofSpace( elements )

cons = dofs.createConstrainer()

cd = dofs.getForTypes( [1,2,3],['u','v'] )

for ido in cd:
  cons.addConstraint( ido , 0.0 ,"main")
  
cons.flush()


###################################
# Store in global data dictionary #
###################################

from pyfem.util.dataStructures import GlobalData

globdat = GlobalData( nodes, elements, dofs )

################################
# Solution procedure (Box 2.3) #
################################

from numpy import zeros, array
from pyfem.fem.Assembly import assembleTangentStiffness

#################################
# Step 1:                       #
# Initialize (Delta a)          #
#################################

a    = globdat.state
Da   = globdat.Dstate
fint = zeros( len(dofs) ) 
fext = zeros( len(dofs) ) 

loadDof = dofs.getForType(4,'v')
Dfext   = zeros( len(dofs) )
Dfext[loadDof] = -2.*DF

output = [ [0.,0.] ]

#Load step iterator
for i in range(N):
  
  print('=================================')
  print(' Load step %i' % i)
  print('=================================')
  print('  NR iter : L2-norm residual')

  #############################################  
  # Step 2:                                   #
  # Compute the new external force vector     #
  #############################################

  fext = fext + Dfext
  
  #Initialize Newton-Raphson iteration parameters  
  error = 1.
  iiter = 0

  #############################################
  # Step 3:                                   #
  # Compute the tangent stiffness matrix      #
  #############################################

  K,fint = assembleTangentStiffness( props, globdat )

  while error > tol:
    
    ###############################################  
    # Step 4+5:                                   #  
    # Solve for da (while satisfying constraints) #
    ###############################################

    da = dofs.solve( K, fext-fint )
    
    ###############################################
    # Step 6:                                     #
    # Update Delta a                              #
    ###############################################

    Da[:] += da[:]
    a [:] += da[:]

    ###############################################
    # Step 7-10                                   #
    # -> Compute Delta epsilon                    #
    # -> Compute Delta sigma                      #
    # -> Compute new sigma                        #
    # -> Compute new internal force vector        #
    ###############################################

    K,fint = assembleTangentStiffness( props, globdat )
    
    ###############################################
    # Step 11:                                    #
    # Convergence check                           #
    ###############################################

    error  = dofs.norm( fext-fint )               

    #Increment the Newton-Raphson iteration counter
    iiter += 1

    print('  Iter', iiter, ':', error)

    if iiter == iterMax:
      raise RuntimeError('Newton-Raphson iterations did not converge!')

  #Update the state vector

  Da[:]  = zeros( len(dofs) )

  #Commit history values
  elements.commitHistory()

  #Store the output
  output.append( [ a[loadDof], fint[loadDof] ] )

  print('=================================')
  

###############################
# Post-processing             #
###############################

from pylab import plot, show, xlabel, ylabel

plot( [-x[0] for x in output], [-0.5*x[1] for x in output], 'ro' )

#Exact solution
from math import sqrt
from numpy import arange

l = lambda v : sqrt( b**2 +(h-v)**2 )
F = lambda v : -E * Area * (h-v)/l(v) * (l(v)-l(0))/l(0) + k * v

vrange = arange(0,1.2,0.01)
plot( vrange, [F(vval) for vval in vrange], 'b-' ) 
xlabel('v [m]')
ylabel('F [N]')

show()
