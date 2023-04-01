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

from numpy import zeros, array, pi
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness

from pyfem.util.logger   import getLogger

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class BuckEigSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol     = 1.0e-3
    self.iterMax = 10 

    BaseModule.__init__( self , props )

    self.fext  = zeros( len(globdat.dofs) )  
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):
      
    K0,fint  = assembleTangentStiffness( props, globdat )

    globdat.state = globdat.dofs.solve( K0 , globdat.fhat )
         
    K,fint  = assembleTangentStiffness( props, globdat )

    eigenvals , eigenvecs = globdat.dofs.eigensolve( K0 , (K0-K) )

    globdat.state = eigenvecs
  
    globdat.elements.commitHistory()

    globdat.active = False 

    self.printResults( eigenvals )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printResults( self , eigenvals):

    logger.info("Staggered solver ............")
    logger.info("    =============================================")
    logger.info("    Eigen modes")
    logger.info("    =============================================")
    logger.info("    Mode   Load")
    
    for i,f in enumerate(eigenvals):
      logger.info(' %3i : %6.4e  ' %(i+1,f))
      
    logger.info('======================================\n')
