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

from numpy import zeros, array, dot
from pyfem.fem.Assembly import assembleTangentStiffness

from pyfem.util.logger   import getLogger

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class RiksSolver( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol       = 1.0e-5
    self.optiter   = 5
    self.iterMax   = 10
    self.fixedStep = False
    self.maxFactor = 1.0e20

    globdat.totalFactor = 1.0
    globdat.factor    = 1.0
    self.maxLam    = 1.0e20

    dofCount    = len(globdat.dofs)

    BaseModule.__init__( self , props )

    if not hasattr(globdat,"Daprev"):
      globdat.Daprev    = zeros( dofCount )
      globdat.Dlamprev  = 1.0

    globdat.lam    = 1.0

    print("\n  Starting Riks arclength solver\n")

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def run( self , props , globdat ):

    stat = globdat.solverStatus
    
    stat.increaseStep()
   
    a    = globdat.state
    Da   = globdat.Dstate
    fhat = globdat.fhat
 
    self.printHeader( stat.cycle )
      
    # Initialize Newton-Raphson iteration parameters  

    error = 1.
 
    # Predictor

    if stat.cycle == 1:    
      K,fint = assembleTangentStiffness( props, globdat )      
      Da1    = globdat.dofs.solve( K , globdat.lam*fhat )
      Dlam1  = globdat.lam
    else:
#      Da1    = self.factor * self.Daprev
#      Dlam1  = self.factor * self.Dlamprev

      Da1    = globdat.factor * globdat.Daprev
      Dlam1  = globdat.factor * globdat.Dlamprev
      globdat.lam += Dlam1
  
    a [:] += Da1[:]
    Da[:] =  Da1[:]

    Dlam = Dlam1

    K,fint = assembleTangentStiffness( props, globdat )  

    res = globdat.lam*fhat-fint  

    while error > self.tol:

      stat.iiter += 1

      d1 = globdat.dofs.solve( K , fhat )
      d2 = globdat.dofs.solve( K , res )
       
      ddlam = -dot(Da1,d2)/dot(Da1,d1)
      dda   = ddlam*d1 + d2
       
      Dlam        += ddlam
      globdat.lam += ddlam
      
      Da[:] += dda[:]
      a [:] += dda[:]

      K,fint = assembleTangentStiffness( props, globdat )

      res = globdat.lam*fhat-fint
 
      error  = globdat.dofs.norm( res ) / globdat.dofs.norm( globdat.lam*fhat )

      self.printIteration( stat.iiter,error)

      if stat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')

    # Converged

    self.printConverged( stat.iiter )
    
    globdat.elements.commitHistory()

    globdat.fint = fint

    if not self.fixedStep:
      globdat.factor = pow(0.5,0.25*(stat.iiter-self.optiter))
      globdat.totalFactor *= globdat.factor

    if globdat.totalFactor > self.maxFactor:
      globdat.factor = 1.0

    globdat.Daprev[:] = Da[:]
    globdat.Dlamprev  = Dlam

    if globdat.lam > self.maxLam or stat.cycle > 1000:
      globdat.active=False

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printHeader( self , cycle):
    
    logger.info("Riks solver .................")
    logger.info("    =============================================")
    logger.info("    Load step %i"%cycle)
    logger.info("    =============================================")
    logger.info('    Newton-Raphson   : L2-norm residual')    

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printIteration( self , iiter , error ):

    logger.info('    Iteration %4i   : %6.4e'%(iiter,error) )

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printConverged( self , iiter ):

    logger.info('    ---------------------------------------------')
    logger.info('    Converged in %i iterations' %iiter)

