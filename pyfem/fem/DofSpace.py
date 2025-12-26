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

import numpy as np
from numpy import array,  where
import scipy.linalg

from scipy.sparse.linalg import spsolve, eigsh
from pyfem.util.itemList import itemList
from pyfem.util.fileParser import readNodeTable
from pyfem.util.logger import getLogger
from pyfem.fem.Constrainer import Constrainer

from copy import deepcopy
from typing import Any, List, Sequence, Optional, Tuple

logger = getLogger()


class DofSpace:
    """Representation of the global degrees-of-freedom space.

    The class maps node identifiers and DOF types to global DOF indices and
    provides utility routines for constraint handling, solves and eigenvalue
    computations in the constrained subspace.
    """

    def __init__(self, elements: Any) -> None:
        """Create a DofSpace from an elements container.

        Args:
            elements: Element container providing ``nodes`` and ``getDofTypes()``.
        """

        self.dofTypes = elements.getDofTypes()
        self.dofs = array(list(range(len(elements.nodes) * len(self.dofTypes)))).reshape(
            (len(elements.nodes), len(self.dofTypes))
        )
        self.nodes = elements.nodes

        # Create the ID map
        self.IDmap = itemList()
        for ind, ID in enumerate(elements.nodes):
            self.IDmap.add(ID, ind)

        self.allConstrainedDofs: List[int] = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __str__(self) -> str:
        """Return a string representation of the DOF array."""

        return str(self.dofs)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __len__(self) -> int:
        """Return the total number of global DOFs."""

        return len(self.dofs.flatten())

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def setConstrainFactor(self, fac: float, loadCase: str = "All_") -> None:
        """Set constraint scaling factor for all or a specific load case."""

        if loadCase == "All_":
            for name in self.cons.constrainedFac.keys():
                self.cons.constrainedFac[name] = fac
        else:
            self.cons.constrainedFac[loadCase] = fac

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def readFromFile(self, fname: str) -> None:
        """Read constraint definitions from a file and create a Constrainer.

        The file is parsed with :func:`pyfem.util.fileParser.readNodeTable` and
        passed to :meth:`createConstrainer`.
        """

        logger.info("  Reading constraints")
        logger.info("  -----------------------------------------------------------")

        nodeTable = readNodeTable(fname, "NodeConstraints", self.nodes)

        self.cons = self.createConstrainer(nodeTable)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
              
    def createConstrainer(self, nodeTables: Optional[Sequence[Any]] = None) -> Constrainer:
        """Create and return a :class:`Constrainer` from parsed node tables.

        If ``nodeTables`` is None a default main constraint group is created.
        Otherwise the function iterates over provided node tables and registers
        constraints accordingly.
        """

        cons = Constrainer(len(self))

        if nodeTables is None:

            label = "main"
            cons.constrainedDofs[label] = []
            cons.constrainedVals[label] = []
            cons.constrainedFac[label] = 1.0

            self.cons = cons
            return cons

        for nodeTable in nodeTables:

            label = nodeTable.subLabel

            cons.constrainedDofs[label] = []
            cons.constrainedVals[label] = []
            cons.constrainedFac[label] = 1.0

            for item in nodeTable.data:

                nodeID = item[1]
                dofType = item[0]
                val = item[2]

                if not nodeID in self.nodes:
                    raise RuntimeError("Node ID " + str(nodeID) + " does not exist")

                ind = self.IDmap.get(nodeID)

                if dofType not in self.dofTypes:
                    raise RuntimeError('DOF type "' + dofType + '" does not exist')

                if len(item) == 3:
                    dofID = self.dofs[ind, self.dofTypes.index(dofType)]

                    cons.addConstraint(dofID, val, label)
                else:
                    slaveNodeID = item[4]
                    slaveDofType = item[3]
                    factor = item[5]

                    if not slaveNodeID[0] in self.nodes:
                        raise RuntimeError("Node ID " + str(slaveNodeID) + " does not exist")

                    slaveInd = self.IDmap.get(slaveNodeID)

                    if slaveDofType not in self.dofTypes:
                        raise RuntimeError('DOF type "' + slaveDofType + '" does not exist')

                    slaveDof = self.dofs[slaveInd, self.dofTypes.index(slaveDofType)]

                    dofID = self.dofs[ind, self.dofTypes.index(dofType)]

                    cons.addConstraint(dofID, [val, slaveDof, factor], label)

        # Check for all tyings whether master of slave is not slave itself
        cons.checkConstraints(self, nodeTables)

        cons.flush()

        return cons

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getForType(self, nodeIDs: Sequence[Any], dofType: str) -> array:
        """Return DOF indices for a given DOF type and list of node IDs."""

        return self.dofs[self.IDmap.get(nodeIDs), self.dofTypes.index(dofType)]
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getForTypes(self, nodeIDs: Sequence[Any], dofTypes: Sequence[str]) -> List[int]:
        """Return DOF indices for multiple DOF types and multiple nodes."""

        dofs: List[int] = []

        for node in nodeIDs:
            for dofType in dofTypes:
                dofs.append(self.dofs[self.IDmap.get(node), self.dofTypes.index(dofType)])

        return dofs
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getDofName(self, dofID: int) -> str:
        """Return a human readable name for a DOF, e.g. 'u[14]'."""

        return self.getTypeName(dofID) + '[' + str(self.getNodeID(dofID)) + ']'

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
    def getNodeID(self, dofID: int) -> Any:
        """Return the node identifier associated with a DOF index."""

        return self.nodes.findID(int(where(self.dofs == dofID)[0]))
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getType(self, dofID: int) -> int:
        """Return the local DOF type index for a global DOF id."""

        return int(where(self.dofs == dofID)[1])
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getTypeName(self, dofID: int) -> str:
        """Return the DOF type name for a given DOF id."""

        return self.dofTypes[self.getType(dofID)]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def get(self, nodeIDs: Sequence[Any]) -> array:
        """Return all DOF ids for the provided node IDs."""

        return self.dofs[self.IDmap.get(nodeIDs)].flatten()

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def copyConstrainer(self, dofTypes: Optional[Sequence[str]] = None) -> Constrainer:
        """Return a copy of the current constrainer with additional DOF types."""

        newCons = deepcopy(self.cons)

        if type(dofTypes) is str:
            dofTypes = [dofTypes]

        if dofTypes is None:
            return newCons

        for dofType in dofTypes:
            for iDof in self.dofs[:, self.dofTypes.index(dofType)]:
                for label in newCons.constrainedFac.keys():
                    newCons.addConstraint(iDof, 0.0, label)

        newCons.flush()

        return newCons

#-------------------------------------------------------------------------------
#  
#-------------------------------------------------------------------------------

    def solve(self, A: array, b: array, constrainer: Optional[Constrainer] = None) -> array:
        """Solve the linear system Ax=b respecting constraints and return x.

        For a matrix problem the constrained system is assembled and solved
        in the reduced space. For a diagonal 'A' (len(A.shape)==1) the solve
        is performed element-wise.
        """

        if constrainer is None:
            constrainer = self.cons

        if len(A.shape) == 2:

            a = np.zeros(len(self))

            constrainer.addConstrainedValues(a)

            A_constrained = constrainer.C.transpose() * (A * constrainer.C)

            b_constrained = constrainer.C.transpose() * (b - A * a)

            x_constrained = spsolve(A_constrained, b_constrained)

            x = constrainer.C * x_constrained

            constrainer.addConstrainedValues(x)

        elif len(A.shape) == 1:
            x = b / A

            constrainer.setConstrainedValues(x)

        return x
    
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def eigensolve(self, A: array, B: array, count: int = 5) -> Tuple[array, array]:
        """Compute the lowest ``count`` eigenpairs for the generalized problem A x = lambda B x.

        The computation is performed in the constrained subspace and the eigenvectors
        are expanded back to the full DOF space before returning.
        """

        A_constrained = self.cons.C.T @ ( A @ self.cons.C )
        B_constrained = self.cons.C.T @ ( B @ self.cons.C )

        eigvals, eigvecs = eigsh(A_constrained, count, B_constrained, sigma=0.0, which="LM")

        x = np.zeros(shape=(len(self), count))

        for i, psi in enumerate(eigvecs.transpose()):
            x[:, i] = self.cons.C * psi

        return eigvals, x

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def norm(self, r: array, constrainer: Optional[Constrainer] = None) -> float:
        """Return the norm of ``r`` excluding constrained DOFs."""

        if constrainer is None:
            constrainer = self.cons

        return scipy.linalg.norm(constrainer.C.transpose() * r)
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def maskPrescribed(self, a: array, val: float = 0.0, constrainer: Optional[Constrainer] = None) -> array:
        """Replace prescribed DOFs in ``a`` with ``val`` and return the array."""

        if constrainer is None:
            constrainer = self.cons

        a[constrainer.constrainedDofs["None"]] = val

        return a
