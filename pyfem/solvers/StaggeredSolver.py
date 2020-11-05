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
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness, commit
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
    self.dtime    = 0.1
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
          
    logger.info("Starting staggered solver .......")
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):

    stat = globdat.solverStatus
    stat.increaseStep()

    fext  = zeros( len(globdat.dofs) ) 
    
    for solver in self.solvers:
           
      globdat.iiter = 0
      error         = 1.0
      
      K,fint = assembleTangentStiffness( props, globdat )
       
      solver.cons.setConstrainFactor(1.0)
      
      da = globdat.dofs.solve( K, fext - fint, solver.cons )
      
      globdat.state += da  
     
      if solver.type == "Nonlinear":
      
        norm = globdat.dofs.norm( fext - fint )
      
        while error > self.tol:
        
          stat.iiter += 1
          
          K,fint = assembleTangentStiffness( props, globdat )
       
          solver.cons.setConstrainFactor(0.0)
          
          da = globdat.dofs.solve( K, fext - fint, solver.cons )

          globdat.state += da  
          
          if norm < 1.0e-16:
            error = globdat.dofs.norm( fext-fint )
          else:
            error = globdat.dofs.norm( fext-fint ) / norm
      
          if globdat.iiter == self.iterMax:
            raise RuntimeError('Newton-Raphson iterations did not converge!')
         
    # Combine results and calculate stresses
    
    globdat.fint = assembleInternalForce( props, globdat )

    commit ( props, globdat )    

    globdat.elements.commitHistory()
    
    if stat.cycle == self.maxCycle: # or globdat.lam > self.maxLam:
      globdat.active = False 
    
