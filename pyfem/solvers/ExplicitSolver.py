# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util.BaseModule import BaseModule
from time import time

from numpy import zeros, array, dot
from pyfem.fem.Assembly import assembleInternalForce, assembleMassMatrix

import sys

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class ExplicitSolver ( BaseModule ):

  def __init__( self , props , globdat ):
    
    self.maxCycle = sys.maxsize

    BaseModule.__init__( self , props )
    
    M,self.Mlumped = assembleMassMatrix( props , globdat )
 
    self.loadfunc = eval ( "lambda t : " + str(self.lam) )
    
    globdat.solverStatus.dtime = self.dtime


#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------


  def run( self , props , globdat ):
    
    stat = globdat.solverStatus
    
    stat.increaseStep()
 
    lam  = self.loadfunc( stat.time )
    
    disp = globdat.state
    velo = globdat.velo
    acce = globdat.acce

    fint = globdat.fint
    fhat = globdat.fhat
    
    velo += 0.5*stat.dtime * acce;
    disp += stat.dtime * velo
    
    fint  = assembleInternalForce( props, globdat )

    globdat.dofs.setConstrainFactor(lam)

    acce = globdat.dofs.solve( self.Mlumped , lam*fhat-fint )
       
    velo += 0.5 * stat.dtime * acce

    globdat.acce[:] = acce[:]
  
    globdat.elements.commitHistory()

    self.printStep( globdat )
    
    if stat.cycle == self.maxCycle:
      globdat.active = False

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

  def printStep( self , globdat ):
 
    stat = globdat.solverStatus
    
    if stat.cycle%20 == 0 or stat.cycle == 1:
      print("  Cycle     Time         Kin.Energy")
      print("  ---------------------------------------")
  
    print(f' {stat.cycle:5d} ', end=' ')
    print(f' {stat.time:10.3e} ', end=' ')
    print(f' {0.5 * dot(globdat.velo, (self.Mlumped * globdat.velo)):10.3e} ')

