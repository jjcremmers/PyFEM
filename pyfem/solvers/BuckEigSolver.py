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

from pyfem.util.BaseModule import BaseModule
from pyfem.util.dataStructures import Properties, GlobalData

from numpy import zeros, array, pi
from pyfem.fem.Assembly import assembleInternalForce, assembleTangentStiffness
from pyfem.fem.Assembly import assembleExternalForce

from pyfem.util.logger   import getLogger

from typing import Any, Sequence, Tuple

logger = getLogger()

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

class BuckEigSolver(BaseModule):
    """
    Eigenvalue solver for buckling analysis.

    This solver assembles the initial and current tangent stiffness matrices,
    solves the static problem to obtain the prebuckling state and computes
    eigenvalues/eigenvectors for the buckling problem.

    Attributes:
        tol: convergence tolerance for internal iterations
        iterMax: maximum iterations (unused currently)
        fext: external force vector (numpy array-like)
    """

    def __init__(self, props: Properties, globdat: GlobalData) -> None:
        """
        Initialize solver.

        Args:
            props: analysis properties
            globdat: global data structure containing mesh, dofs and state
        """
        self.tol     = 1.0e-3
        self.iterMax = 10

        BaseModule.__init__(self, props)

        # external force vector sized to global dof vector
        self.fext: Any = zeros(len(globdat.dofs))

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def run(self, props: Properties, globdat: GlobalData) -> None:
        """
        Execute the buckling eigenvalue analysis.

        Steps:
        - assemble initial tangent stiffness K0 and internal force fint
        - assemble external force fext and solve for prebuckling state
        - assemble updated stiffness K and compute eigenpairs using DOF solver
        - commit element history and deactivate model if appropriate

        Args:
            props: analysis properties
            globdat: global data structure (modified in-place)
        """
        self.writeHeader()

        # initial tangent stiffness and internal force
        K0, fint = assembleTangentStiffness(props, globdat)  # type: Tuple[Any, Any]

        # assemble external force
        fext = assembleExternalForce(props, globdat)  # type: Any

        # solve static problem for prebuckling state
        globdat.state = globdat.dofs.solve(K0, fext)

        # assemble updated stiffness at prebuckling state
        K, fint = assembleTangentStiffness(props, globdat)  # type: Tuple[Any, Any]

        # compute eigenvalues/eigenvectors for buckling: generalized problem
        globdat.eigenvals, globdat.eigenvecs = globdat.dofs.eigensolve(K0, (K0 - K))

        # commit element history and finalize
        globdat.elements.commitHistory()
        globdat.active = False

        self.printResults(globdat.eigenvals)

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

    def printResults(self, eigenvals: Sequence[float]) -> None:
        """
        Print computed eigenvalues.

        Args:
            eigenvals: iterable of eigenvalues (loads) to report
        """
        logger.info("   Eigen modes")
        logger.info("   ---------------------------------------------------------")
        logger.info("    Mode |  Load")

        for i, f in enumerate(eigenvals):
            logger.info(f"    {i+1:4d} | {f:6.4e}  ")
