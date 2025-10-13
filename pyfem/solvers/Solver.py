################################################################################
#  This Python file is part of PyFEM, the code that accompanies the book:      #
#                                                                              #
#    'Non-Linear Finite Element Analysis of Solids and Structures'             #
#    R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel            #
#    John Wiley and Sons, 2012, ISBN 978-0470666449                            #
#                                                                              #
#  Copyright (C) 2011-2025. The code is written in 2011-2012 by                #
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

class Solver:
    def __init__(self, props, globdat):
        solverProps = getattr(props, "solver")
        solverType  = solverProps.type  # e.g. "DissipatedEnergySolver"

        try:
            module = import_module(f"pyfem.solvers.{solverType}")
        except ModuleNotFoundError as e:
            raise ImportError(
                f"Solver module 'pyfem.solvers.{solverType}' not found. "
                f"Check the 'type' in your input file."
            ) from e

        try:
            solver_cls = getattr(module, solverType)
        except AttributeError as e:
            raise ImportError(
                f"Class '{solverType}' not found in module 'pyfem.solvers.{solverType}'. "
                f"Ensure the class name matches the file name."
            ) from e

        props.currentModule = "solver"
        self.solver = solver_cls(props, globdat)
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
        
  def run( self , props , globdat ):

    self.solver.run( props , globdat )
