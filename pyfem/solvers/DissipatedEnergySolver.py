# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

################################################################################
#  This solver is presented in detail in the following paper:                  #
#                                                                              #
#  E.Borjesson, J.J.C. Remmers, M. Fagerstrom (2023) A generalised             #  
#    path-following solver for robust analysis of material failure,            #
#    Computational Mechanics, doi: 10.1007/s00466-022-02175-w                  #
#                                                                              #
################################################################################

from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array, dot
from pyfem.fem.Assembly import assembleTangentStiffness,assembleDissipation,assembleExternalForce

from pyfem.util.logger   import getLogger

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class DissipatedEnergySolver( BaseModule ):

  def __init__( self , props , globdat ):

    self.tol       = 1.0e-4
    self.optiter   = 5
    self.iterMax   = 10
    self.maxdTau   = 1.0e20

    self.factor    = 1.0
    self.maxLam    = 1.0e20
    self.lam       = 1.0
    self.disstype  = "Local"

    dofCount    = len(globdat.dofs)

    BaseModule.__init__( self , props )

    self.method    = "arclength-controlled"
    self.Dlam      = self.lam

    globdat.lam    = self.lam
    globdat.dTau   = 0.0

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def run( self , props , globdat ):

    stat = globdat.solverStatus    
    stat.increaseStep()
    
    self.writeHeader( stat.cycle )
   
    a    = globdat.state
    Da   = globdat.Dstate
    
    fhat = assembleExternalForce( props, globdat )  
 
    self.printHeader( stat.cycle )
      
    error         = 1.
    lam0          = globdat.lam
    
    if stat.cycle == 1:       
      K,fint = assembleTangentStiffness( props, globdat )      
      Da1    = globdat.dofs.solve( K , globdat.lam*fhat )
      Dlam   = globdat.lam
    else:
      Da1    = self.factor * self.Daprev
      Dlam   = self.factor * self.Dlamprev
      globdat.lam += Dlam

    a [:] += Da1[:]
    Da[:] =  Da1[:]
            
    K,fint    = assembleTangentStiffness( props, globdat )  
    dgda,diss = assembleDissipation( props , globdat )   

    res = globdat.lam*fhat-fint  

    while error > self.tol:
      stat.iiter += 1

      d1 = globdat.dofs.solve( K , fhat )
      d2 = globdat.dofs.solve( K , res )
            
      ddlamR = -dot(Da1,d2)/dot(Da1,d1)
      ddaR   = ddlamR*d1 + d2

      if self.method == 'nrg-controlled':
        if self.disstype == "Classic":
          h  =  0.5 * lam0 * fhat  # wordt dgda
          w  = -0.5 * dot ( (a-Da) , fhat )   # wordt 0.0
          g  =  0.5 * dot ( ( lam0 * Da - Dlam * ( a[:] - Da[:] ) ) , fhat ) - globdat.dtau  
        elif self.disstype == "Local":
          h = dgda
          w = 0.0
          g = diss - globdat.dtau
        else:
          raise RuntimeError('Not implemented!')
  
        denom  = - dot ( h , d1 ) - w

        ddaN    = d2 - ( -d1 * ( dot( h , d2 ) + g ) ) / denom
        ddlamN  = -g - ( dot( -1.0*h , d2 ) - g * ( 1.0 + denom ) ) / denom;
 
      if self.method == 'arclength-controlled':
        Dlam        += ddlamR
        globdat.lam += ddlamR
            
        Da[:] += ddaR[:]
        a [:] += ddaR[:]        
      elif self.method == 'nrg-controlled':
        Dlam        += ddlamN
        globdat.lam += ddlamN
            
        Da[:] += ddaN[:]
        a [:] += ddaN[:]
                              
      # Solve for new displacement vector, load factor      

      dgda,diss = assembleDissipation( props , globdat )         
      K,fint = assembleTangentStiffness( props, globdat )
    
      res = globdat.lam*fhat-fint
      # Check convergence

      error  = globdat.dofs.norm( res ) / globdat.dofs.norm( globdat.lam*fhat )

      self.printIteration( stat.iiter , error )

    # If converged, calculate the amount of energy that has been dissipated in the \
    # previous step.

    if self.disstype == "Classic":
      diss = 0.5 * dot( ( lam0 * Da - Dlam * ( a[:] - Da[:] ) ),fhat )
   

    self.printConverged( stat.iiter , diss )

    self.Daprev    = Da
    self.Dlamprev  = Dlam
  
    self.factor = pow(0.5,0.25*(stat.iiter-self.optiter))
  
    globdat.dtau = diss
        
    if self.method == 'arclength-controlled':
      if diss > self.switchEnergy:
        print('   Switch to nrg diss. arc-length')
        self.method       = 'nrg-controlled'
        globdat.dtau = diss
        self.factor  = 1.0     
    else:
      globdat.dtau *= self.factor
      
      if globdat.dtau > self.maxdTau:
        self.factor *= self.maxdTau / globdat.dtau
        globdat.dtau = self.maxdTau
            
    globdat.elements.commitHistory()

    globdat.fint = fint

    if globdat.lam > self.maxLam or stat.cycle > self.maxCycle:
      globdat.active=False
      
    self.writeFooter( globdat )       
 
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printHeader( self , cycle):

    logger.info('    Newton-Raphson   : L2-norm residual')

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printIteration( self , iiter , error ):

    logger.info(f'    Iteration {iiter:4d}   : {error:6.4e}')

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printConverged( self , iiter , dissnrg ):

    logger.info('    ---------------------------------------------')
    logger.info(f'    Converged in {iiter} iterations with {self.method}.\n')
