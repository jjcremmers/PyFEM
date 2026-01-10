# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

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

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def run( self , props , globdat ):
    
    stat = globdat.solverStatus    
    stat.increaseStep()
    
    self.writeHeader( stat.cycle )     
   
    a    = globdat.state
    Da   = globdat.Dstate
    fhat = globdat.fhat
 
    logger.info('    Newton-Raphson............ : L2-norm residual') 
      
    # Initialize Newton-Raphson iteration parameters  

    error = 1.
 
    # Predictor

    if stat.cycle == 1:    
      K,fint = assembleTangentStiffness( props, globdat )      
      Da1    = globdat.dofs.solve( K , globdat.lam*fhat )
      Dlam1  = globdat.lam
    else:
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

      logger.info(f'    Iteration {stat.iiter:4d} ........... : {error:6.4e}')

      if stat.iiter == self.iterMax:
        raise RuntimeError('Newton-Raphson iterations did not converge!')

    # Converged

    logger.info('                                 Converged')
    
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
      
    self.writeFooter( globdat )      
