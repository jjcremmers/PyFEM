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

    print("\n  Starting explicit solver .....\n")

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
  
    print(' %5i ' % stat.cycle, end=' ')
    print(' %10.3e ' % stat.time, end=' ')  
    print(' %10.3e ' % float(0.5*dot(globdat.velo,(self.Mlumped*globdat.velo))))
