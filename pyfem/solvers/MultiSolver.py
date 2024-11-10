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

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

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
              
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def run( self , props , globdat ):
            
        self.solverList[self.iSlv].run( props, globdat )
        
        if not globdat.active and self.iSlv < len(self.solverList)-1:
            self.iSlv += 1
            globdat.active = True
