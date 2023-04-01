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
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness, commit
from pyfem.fem.Assembly import assembleExternalForce
from pyfem.util.logger import getLogger
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

    self.fext  = zeros( len(globdat.dofs) )  
    
    self.solvers = []
        
    dofsToRemove = [e for e in globdat.dofs.dofTypes if e not in self.solver1.dofTypes ]
    self.solver1.cons = globdat.dofs.copyConstrainer(dofsToRemove)     
    
    dofsToRemove = [e for e in globdat.dofs.dofTypes if e not in self.solver2.dofTypes ]        
    self.solver2.cons = globdat.dofs.copyConstrainer(dofsToRemove)
    
    self.solvers.append(self.solver1)
    self.solvers.append(self.solver2)
    
    self.loadfunc = eval ( "lambda t : " + str(self.loadFunc) )
    globdat.solverStatus.dtime = self.dtime
           
    logger.info("Starting staggered solver .......")
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):

    stat = globdat.solverStatus
    
    self.stat = stat
    
    stat.increaseStep()

    #fext  = zeros( len(globdat.dofs) ) 
    
    logger.info("Staggered solver ............")
    logger.info("    =============================================")
    logger.info("    Load step %i"%globdat.solverStatus.cycle)
    logger.info("    =============================================")
    
    for solver in self.solvers:
           
      stat.iiter = 0
      error      = 1.0
      
      K,fint = assembleTangentStiffness( props, globdat )
      fext   = assembleExternalForce   ( props, globdat )       
      
      self.setLoadAndConstraints( solver.cons )
      
      da = globdat.dofs.solve( K, fext - fint, solver.cons )
      
      globdat.state += da  
      
      logger.info('    Solver           : %s' %solver.name )
          
      if solver.type == "Nonlinear":
      
        if solver.name == "dummy":
          norm = 1.0
        else:
          norm = globdat.dofs.norm( fext - fint, solver.cons )
        
        logger.info('    Newton-Raphson   : L2-norm residual')
      
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
            
          logger.info('    Iteration %4i   : %6.4e'%(stat.iiter,error) )            
      
          if stat.iiter == self.iterMax:
            raise RuntimeError('Newton-Raphson iterations did not converge!')
         
    # Combine results and calculate stresses
    
    globdat.fint = assembleInternalForce( props, globdat )

    commit ( props, globdat )    

    globdat.elements.commitHistory()
    
    if stat.cycle == self.maxCycle: # or globdat.lam > self.maxLam:
      globdat.active = False 

#---------------------------------------------------------------------------
#
#  -------------------------------------------------------------------------


  def setLoadAndConstraints( self , cons ):

    #logger.info("    Load step %i"%self.stat.cycle)
 
    time  = self.stat.time
    time0 = time - self.stat.dtime
    
    lam     = self.loadfunc( time  )
    lam0    = self.loadfunc( time0 )

    dlam = lam - lam0
    
    cons.setConstrainFactor( dlam )

    #logger.info('  ---- main load --------------------\n-----')
    #logger.info('    loadFactor       : %4.2f'%lam)
    #logger.info('    incr. loadFactor : %4.2f'%dlam)

    for loadCase in self.loadCases:
      loadProps = getattr( self.myProps, loadCase )
      
      loadfunc = eval ( "lambda t : " + str(loadProps.loadFunc) )
      lam  = loadfunc( time  )
      lam0 = loadfunc( time0 )
      dlam = lam - lam0
      cons.setConstrainFactor( dlam , loadProps.nodeTable )

      #print('  ---- ',loadCase,' ---------------------')
      #print('    loadFactor       : %4.2f'%lam)
      #logger.info('    incr. loadFactor : %4.2f'%dlam)   
