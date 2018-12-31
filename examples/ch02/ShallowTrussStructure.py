#############################
# Define problem parameters #
# see Figure 1.1            #
#############################

#Problem dimensions
b = 10. #m
h = 0.5 #m

#Material parameters
k   = 1000. #N/m
EA0 = 5e6   #N/m^2

#External load increment
DF = 50 #N

#Solver parameters
N       = 30
tol     = 1e-6
iterMax = 5

##############################
# Store data in properties   #
##############################

from pyfem.utils.dataStructures import Properties

props = Properties()
props.Trusses = Properties( { 'EA0' : EA0 } )
props.Spring  = Properties( { 'k' : 2*k } )

#############################
# Defining finite element   #
# data structures           #
#############################

#NodeSet
from pyfem.fem.NodeSet import NodeSet

nodes = NodeSet()
nodes.add( 1, [0.,0.] ) #Attachement point of spring (can be positioned anywhere along the y-axis)
nodes.add( 2, [-b,0.] ) #Left support
nodes.add( 3, [ b,0.] ) #Right support
nodes.add( 4, [0.,h ] ) #Loading point

#ElementSet
from pyfem.fem.ElementSet import ElementSet
from pyfem.elements.Truss import Truss
from pyfem.elements.Spring import Spring

elements = ElementSet( nodes )
elements.add( 1, Truss ( [2,4] ) )
elements.add( 2, Truss ( [3,4] ) )
elements.add( 3, Spring( [1,4] ) )

#Add groups
elements.addGroup( 'Trusses', [1,2] )
elements.addGroup( 'Spring' , [3]   )

#DofSpace
from pyfem.fem.DofSpace import DofSpace

dofs = DofSpace( elements )

dofs.constrain( 1, ['u','v'] )
dofs.constrain( 2, ['u','v'] )
dofs.constrain( 3, ['u','v'] )

###################################
# Store in global data dictionary #
###################################

from pyfem.utils.dataStructures import GlobalData

globdat = GlobalData( nodes, elements, dofs )

################################
# Solution procedure (Box 2.3) #
################################

from numpy import zeros, array
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness

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

  K = assembleTangentStiffness( props, globdat )

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

    fint = assembleInternalForce( props, globdat )
    K    = assembleTangentStiffness( props, globdat )
    
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
F = lambda v : -EA0 * (h-v)/l(v) * (l(v)-l(0))/l(0) + k * v

vrange = arange(0,1.2,0.01)
plot( vrange, [F(vval) for vval in vrange], 'b-' ) 
xlabel('v [m]')
ylabel('F [N]')

show()
