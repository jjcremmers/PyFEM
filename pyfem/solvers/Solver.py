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

import importlib
from typing import Any

from pyfem.util.dataStructures import Properties, GlobalData


class Solver:
    """
    Loader and base wrapper for solver implementations.

    This class dynamically imports the solver module specified in
    props.solver.type, instantiates the solver class and delegates the
    run call to that instance.

    Attributes:
        solver: Instantiated solver object (implementation-specific)
    """

    def __init__(self, props: Properties, globdat: GlobalData) -> None:
        """
        Load and instantiate the solver specified in properties.

        Args:
            props: Analysis properties (must contain a 'solver' attribute
                   with a 'type' field indicating the solver class name)
            globdat: GlobalData instance passed to the solver constructor

        Raises:
            ImportError: If the solver module or class cannot be found.
        """
        solverProps = getattr(props, "solver")
        solverType: str = solverProps.type

        try:
            mod = importlib.import_module(f"pyfem.solvers.{solverType}")
            SolverClass = getattr(mod, solverType)
        except ModuleNotFoundError as e:
            raise ImportError(
                f"Solver module 'pyfem.solvers.{solverType}' not found."
            ) from e
        except AttributeError as e:
            raise ImportError(
                f"Class '{solverType}' not found in module "
                f"'pyfem.solvers.{solverType}'."
            ) from e

        props.currentModule = "solver"
        self.solver: Any = SolverClass(props, globdat)

    def run(self, props: Properties, globdat: GlobalData) -> None:
        """
        Delegate execution to the selected solver instance.

        Args:
            props: Analysis properties
            globdat: GlobalData instance (may be modified in-place)
        """
        self.solver.run(props, globdat)
