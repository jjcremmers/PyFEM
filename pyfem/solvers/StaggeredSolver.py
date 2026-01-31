# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util import logger
from pyfem.util.BaseModule import BaseModule

import numpy as np
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness, commit
from pyfem.fem.Assembly import assembleExternalForce
from pyfem.util.logger import getLogger,separator
import sys

logger = getLogger()


#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class StaggeredSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol      = 1.0e-3
    self.iterMax  = 10 

    self.maxCycle = sys.maxsize
    self.maxLam   = 1.0e20
    self.dtime    = 1.0
    self.loadFunc = "t"
    self.loadCases= []
    
    BaseModule.__init__( self , props )

    self.fext  = np.zeros( len(globdat.dofs) )  
    
    self.solvers = []
        
    dofsToRemove = [e for e in globdat.dofs.dofTypes if e not in self.solver1.dofTypes ]
    self.solver1.cons = globdat.dofs.copyConstrainer(dofsToRemove)     
    
    dofsToRemove = [e for e in globdat.dofs.dofTypes if e not in self.solver2.dofTypes ]        
    self.solver2.cons = globdat.dofs.copyConstrainer(dofsToRemove)
    
    self.solvers.append(self.solver1)
    self.solvers.append(self.solver2)
    
    self.loadfunc = eval ( "lambda t : " + str(self.loadFunc) )
    globdat.solverStatus.dtime = self.dtime
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):

    stat = globdat.solverStatus
    
    self.stat = stat
    
    stat.increaseStep()
    
    self.writeHeader( stat.cycle )    

    for solver in self.solvers:
           
      stat.iiter = 0
      error      = 1.0
      
      K,fint = assembleTangentStiffness( props, globdat )
      fext   = assembleExternalForce   ( props, globdat )       
      
      self.setLoadAndConstraints( solver.cons )
      
      da = globdat.dofs.solve( K, fext - fint, solver.cons )
      
      globdat.state += da  
      
      logger.info(f'    Solver.................... : {solver.name}')
          
      if solver.type == "Nonlinear":
      
        if solver.name == "dummy":
          norm = 1.0
        else:
          norm = globdat.dofs.norm( fext - fint, solver.cons )
        
        logger.info('      Newton-Raphson........... : L2-norm residual')
      
        while error > self.tol:
        
          stat.iiter += 1
          
          K,fint = assembleTangentStiffness( props, globdat )
       
          solver.cons.setConstrainFactor(0.0)
          
          da = globdat.dofs.solve( K, fext - fint, solver.cons )

          globdat.state += da  
          
          if solver.name == "dummy":
            error = 1.0e-8           
          elif norm < 1.0e16:
            error = globdat.dofs.norm( fext-fint, solver.cons )
          else:
            error = globdat.dofs.norm( fext-fint ) / norm
            
          logger.info(f'      Iteration {stat.iiter:4d} ........... : {error:6.4e}')
      
          if stat.iiter == self.iterMax:
            raise RuntimeError('Newton-Raphson iterations did not converge!')
         
    # Combine results and calculate stresses
    
    globdat.fint = assembleInternalForce( props, globdat )

    commit ( props, globdat )    

    globdat.elements.commitHistory()
    
    if stat.cycle == self.maxCycle: # or globdat.lam > self.maxLam:
      globdat.active = False 

    separator()
    self.writeFooter( globdat )      
      

#---------------------------------------------------------------------------
#
#  -------------------------------------------------------------------------


  def setLoadAndConstraints( self , cons ):
 
    time  = self.stat.time
    time0 = time - self.stat.dtime
    
    lam     = self.loadfunc( time  )
    lam0    = self.loadfunc( time0 )

    dlam = lam - lam0
    
    cons.setConstrainFactor( dlam )

    for loadCase in self.loadCases:
      loadProps = getattr( self.myProps, loadCase )
      
      loadfunc = eval ( "lambda t : " + str(loadProps.loadFunc) )
      lam  = loadfunc( time  )
      lam0 = loadfunc( time0 )
      dlam = lam - lam0
      cons.setConstrainFactor( dlam , loadProps.nodeTable )
