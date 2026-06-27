# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array, pi
from pyfem.fem.Assembly import assembleTangentStiffness, assembleMassMatrix

from pyfem.util.logger   import getLogger
from numpy import sqrt
import h5py

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class DynEigSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol        = 1.0e-3
    self.eigenCount = 5

    BaseModule.__init__( self , props ) 
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
   
  def run( self , props , globdat ):
  
    self.writeHeader()   
      
    K,fint  = assembleTangentStiffness( props, globdat )
         
    M,_ = assembleMassMatrix      ( props , globdat )

    eigenvals , globdat.eigenvecs = globdat.dofs.eigensolve( K , M , self.eigenCount )
    
    globdat.eigenvals = sqrt( eigenvals )
         
    globdat.elements.commitHistory()

    globdat.active = False 
    
    self.printResults( globdat.eigenvals )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printResults( self , eigenvals):

    logger.info('   Eigenfrequencies')
    logger.info("   ----------------------------------------------------------")
    logger.info('   Mode |   Eigenvalue       |  Frequency')
        
    for i,val in enumerate(eigenvals):
      logger.info(f"   {i+1:4d} |   {val:6.4e} rad/s |  {val/(2.0*pi):6.4e} Hz")
