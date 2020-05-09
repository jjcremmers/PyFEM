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

from pyfem.util.logger   import getLogger

logger = getLogger()

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
 
    logger.info("Starting nonlinear solver .........")

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

    logger.info('  NR iter  : L2-norm residual')
    
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

      logger.info('    Iteration %4i   : %6.4e'%(globdat.iiter,error) )

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

    logger.info("    Load step %i"%globdat.cycle)
 
    globdat.lam  = self.loadfunc( globdat.time )
    lam0         = self.loadfunc( globdat.time - self.dtime )

    globdat.dlam = globdat.lam - lam0
    globdat.dofs.setConstrainFactor( globdat.dlam )

    logger.info('  ---- main load --------------------\n-----')
    logger.info('    loadFactor       : %4.2f'%globdat.lam)
    logger.info('    incr. loadFactor : %4.2f'%globdat.dlam)

    for loadCase in self.loadCases:
      loadProps = getattr( self.myProps, loadCase )
      
      loadfunc = eval ( "lambda t : " + str(loadProps.loadFunc) )
      lam  = loadfunc( globdat.time )
      lam0 = loadfunc( globdat.time - self.dtime )
      dlam = lam - lam0
      globdat.dofs.setConstrainFactor( dlam , loadProps.nodeTable )

      logger.info('  ---- ',loadCase,' ---------------------')
      logger.info('    loadFactor       : %4.2f'%lam)
      logger.info('    incr. loadFactor : %4.2f'%dlam)

    fhat  = globdat.fhat
    
    return globdat.lam*fhat

      
