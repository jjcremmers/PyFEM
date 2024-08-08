################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2024. The code is written in 2011-2012 by                #
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
from pyfem.fem.Assembly import assembleTangentStiffness, assembleMassMatrix

from pyfem.util.logger   import getLogger
from math import sqrt
import h5py

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class DynEigSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol        = 1.0e-3
    self.eigenCount = 5
    self.writeToH5  = False

    BaseModule.__init__( self , props ) 
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):
  
    self.writeHeader()   
      
    K,fint  = assembleTangentStiffness( props, globdat )
         
    M,mlump = assembleMassMatrix      ( props , globdat )

    eigenvals , globdat.eigenvecs = globdat.dofs.eigensolve( K , M , self.eigenCount )
    
    globdat.eigenvals = []
    
    for val in eigenvals:
      globdat.eigenvals.append( sqrt(val) )
         
    globdat.elements.commitHistory()

    globdat.active = False 
    
    if self.writeToH5:
      h5file = h5py.File( "modes.h5", 'w')
    
      h5file.create_dataset("modes", globdat.eigenvecs.shape, 
                            dtype='f', data=globdat.eigenvecs)   
  
    self.printResults( globdat.eigenvals )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printResults( self , eigenvals):

    logger.info('   Eigenfrequencies')
    logger.info("   ----------------------------------------------------------")
    logger.info('   Mode |   Eigenvalue       |  Frequency')
        
    for i,val in enumerate(eigenvals):
      logger.info('   %4i |   %6.4e rad/s |  %6.4e Hz' %(i+1,val,val/(2.0*pi)))     
