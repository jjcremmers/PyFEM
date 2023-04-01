################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2023. The code is written in 2011-2012 by                #
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

from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness
from pyfem.fem.Assembly import assembleExternalForce
from math import sin,cos,exp

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
    self.dtime    = 1.0
    self.loadFunc = "t"
    self.loadCases= []

    BaseModule.__init__( self , props )

    if self.maxLam > 1.0e19 and self.maxCycle == sys.maxsize:
      self.maxCycle = 5

    globdat.lam = 0.0
    globdat.solverStatus.dtime = self.dtime

    self.loadfunc = eval ( "lambda t : " + str(self.loadFunc) )
       
    if hasattr(self,"loadTable"):
      self.maxCycle      = len(self.loadTable)
      loadTable          = zeros(self.maxCycle+1)
      loadTable[1:]      = self.loadTable
      self.loadTable     = loadTable
 
    logger.info("Starting nonlinear solver .........")

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def run( self , props , globdat ):

    stat = globdat.solverStatus
    
    stat.increaseStep()
    
    dofCount = len(globdat.dofs)
    
    a     = globdat.state
    Da    = globdat.Dstate

    Da[:] = zeros( dofCount )
    fint  = zeros( dofCount ) 
    
    logger.info("Nonlinear solver ............")
    logger.info("    =============================================")
    logger.info("    Load step %i"%globdat.solverStatus.cycle)
    logger.info("    =============================================")
    logger.info('    Newton-Raphson   : L2-norm residual')
    
    self.setLoadAndConstraints( globdat )
    
    K,fint = assembleTangentStiffness( props, globdat )

    error = 1.

    self.setLoadAndConstraints( globdat )
    
    fext   = assembleExternalForce   ( props, globdat )
        
    while error > self.tol:

      stat.iiter += 1
	      
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

      logger.info('    Iteration %4i   : %6.4e'%(stat.iiter,error) )

      globdat.dofs.setConstrainFactor( 0.0 )

      if stat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')

    # Converged
    
    globdat.elements.commitHistory()

    Da[:]  = zeros( len(globdat.dofs) )

    globdat.fint = fint
    
    if stat.cycle == self.maxCycle or globdat.lam > self.maxLam:
      globdat.active = False 

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def setLoadAndConstraints( self , globdat ):
 
    if hasattr(self,"loadTable"):
      cycle = globdat.solverStatus.cycle
      
      globdat.lam  = self.loadTable[cycle]
      globdat.dlam = self.loadTable[cycle]-self.loadTable[cycle-1]

      globdat.dofs.setConstrainFactor( globdat.dlam )
            
      globdat.solverStatus.lam = globdat.lam
    else:   
      globdat.lam  = self.loadfunc( globdat.solverStatus.time )
      lam0         = self.loadfunc( globdat.solverStatus.time - globdat.solverStatus.dtime )
    
      globdat.dlam = globdat.lam - lam0
      globdat.dofs.setConstrainFactor( globdat.dlam )
    
      globdat.solverStatus.lam = globdat.lam

      logger.debug('  ---- main load -------------------------')
      logger.debug('    loadFactor       : %4.2f'%globdat.lam)
      logger.debug('    incr. loadFactor : %4.2f'%globdat.dlam)

      for loadCase in self.loadCases:
        loadProps = getattr( self.myProps, loadCase )
        loadfunc = eval ( "lambda t : " + str(loadProps.loadFunc) )
        lam  = loadfunc( globdat.solverStatus.time )
        lam0 = loadfunc( globdat.solverStatus.time - globdat.solverStatus.dtime )
        dlam = lam - lam0
        globdat.dofs.setConstrainFactor( dlam , loadProps.nodeTable )
        
        logger.debug('  ---- %s ---------------------' %loadCase)
        logger.debug('    loadFactor       : %4.2f'%lam)
        logger.debug('    incr. loadFactor : %4.2f'%dlam)

      
