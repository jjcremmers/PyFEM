# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util.BaseModule import BaseModule

from numpy import zeros, array
from pyfem.fem.Assembly import assembleInternalForce, assembleExternalForce
from pyfem.fem.Assembly import assembleTangentStiffness, prepare, commit

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class LinearSolver ( BaseModule ):

  def __init__( self , props , globdat ):

    BaseModule.__init__( self , props )

    self.fext  = zeros( len(globdat.dofs) )  
  
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
  def run( self , props , globdat ):
  
    self.writeHeader() 

    globdat.solverStatus.increaseStep()

    prepare( props , globdat )
      
    K,fint = assembleTangentStiffness( props, globdat )
    fext   = assembleExternalForce   ( props, globdat )

    state0 = globdat.state
         
    globdat.state = globdat.dofs.solve( K, fext )
     
    globdat.Dstate = globdat.state - state0

    globdat.fint = assembleInternalForce( props, globdat )

    commit ( props, globdat )    

    globdat.elements.commitHistory()

    globdat.active = False 
    
    self.writeFooter( globdat )
