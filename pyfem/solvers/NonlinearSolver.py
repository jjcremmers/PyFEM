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
from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness
from math import sin

import sys

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class NonlinearSolver( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol      = 1.0e-3
    self.iterMax  = 10 

    self.maxCycle = sys.maxsize
    self.maxLam   = 1.0e20
    self.dtime    = 0.1
    self.loadFunc = "t"
    self.loadCases= []

    BaseModule.__init__( self , props )

    if self.maxLam > 1.0e19 and self.maxCycle == sys.maxsize:
      self.maxCycle = 5

    globdat.lam = 0.0

    self.loadfunc = eval ( "lambda t : " + str(self.loadFunc) )

    print(self.loadCases)
 
    print("\n  Starting nonlinear solver ....\n")

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def run( self , props , globdat ):

    globdat.cycle += 1
    globdat.time  += self.dtime
    
    dofCount = len(globdat.dofs)

    a     = globdat.state
    Da    = globdat.Dstate

    Da[:] = zeros( dofCount )
    fint  = zeros( dofCount ) 

    globdat.iiter = 0 

    K,fint = assembleTangentStiffness( props, globdat )
    
    error = 1.

    fext = self.setLoadAndConstraints( globdat )
    
    while error > self.tol:

      globdat.iiter += 1
	      
      da = globdat.dofs.solve( K, fext - fint )

      Da[:] += da[:]
      a [:] += da[:]

      K,fint = assembleTangentStiffness( props, globdat )
  
      # note that the code is different from the one presented in the book, which
      # is slightly shorter for the sake of clarity.
      # In the case of a prescribed displacement, the external force is zero
      # and hence its norm is zero. In that case, the norm of the residue is not
      # divided by the norm of the external force.

      norm = globdat.dofs.norm( fext )
  
      if norm < 1.0e-16:
        error = globdat.dofs.norm( fext-fint )
      else:
        error = globdat.dofs.norm( fext-fint ) / norm

      print('  Iter', globdat.iiter, ':', error)

      globdat.dofs.setConstrainFactor( 0.0 )

      if globdat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')

    # Converged
    
    globdat.elements.commitHistory()

    Da[:]  = zeros( len(globdat.dofs) )

    globdat.fint = fint
    
    if globdat.cycle == self.maxCycle or globdat.lam > self.maxLam:
      globdat.active = False 

  def setLoadAndConstraints( self , globdat ):

    print('=================================')
    print(' Load step %i' % globdat.cycle)
    print('=================================')

    globdat.lam  = self.loadfunc( globdat.time )
    lam0         = self.loadfunc( globdat.time - self.dtime )

    globdat.dlam = globdat.lam - lam0
    globdat.dofs.setConstrainFactor( globdat.dlam )

    print('---- MAIN LOAD ------')
    print('  loadFactor       : ',globdat.lam)
    print('  incr. loadFactor : ',globdat.dlam)

    for loadCase in self.loadCases:
      loadProps = getattr( self.myProps, loadCase )
      
      loadfunc = eval ( "lambda t : " + str(loadProps.loadFunc) )
      lam  = loadfunc( globdat.time )
      lam0 = loadfunc( globdat.time - self.dtime )
      dlam = lam - lam0
      globdat.dofs.setConstrainFactor( dlam , loadProps.nodeTable )

      print('---- ',loadCase,' ------')
      print('  loadFactor       : ',lam)
      print('  incr. loadFactor : ',dlam)

    fhat  = globdat.fhat
    
    return globdat.lam*fhat



    print('  NR iter : L2-norm residual')

      
