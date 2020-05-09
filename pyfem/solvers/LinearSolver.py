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

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class LinearSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol     = 1.0e-3
    self.iterMax = 10 

    BaseModule.__init__( self , props )

    self.fext  = zeros( len(globdat.dofs) )  
 
    logger.info("Starting linear solver .......")
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):

    globdat.cycle += 1
      
    K,fint = assembleTangentStiffness( props, globdat )

    state0 = globdat.state
         
    globdat.state = globdat.dofs.solve( K, globdat.fhat )
     
    globdat.Dstate = globdat.state - state0

    globdat.fint = assembleInternalForce( props, globdat )

    commit ( props, globdat )    

    globdat.elements.commitHistory()

    globdat.active = False 
