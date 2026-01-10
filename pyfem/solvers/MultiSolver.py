# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from pyfem.util.BaseModule import BaseModule

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class MultiSolver( BaseModule ):

    def __init__( self , props , globdat ):

        self.solverList = []

        BaseModule.__init__( self , props )
    
        for solverName in self.solvers:
            solverProps = getattr( self.myProps , str(solverName) )
            solverType  = solverProps.type
   
            props.currentModule = "solver."+str(solverName)
      
            exec("from pyfem.solvers."+solverType+" import "+solverType)
           
            self.solverList.append( eval(solverType+"( props , globdat )") )
        
        self.iSlv = 0
              
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def run( self , props , globdat ):
            
        self.solverList[self.iSlv].run( props, globdat )
        
        if not globdat.active and self.iSlv < len(self.solverList)-1:
            self.iSlv += 1
            globdat.active = True
