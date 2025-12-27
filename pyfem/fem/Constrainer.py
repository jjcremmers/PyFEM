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

from numpy import array, dot, zeros
import scipy.linalg
from scipy.sparse import coo_matrix
from typing import Any, Dict, List

from pyfem.util.logger import getLogger

logger = getLogger()


class Constrainer:
    """Container for handling multi-point constraints (ties/prescribed dofs).

    The class records constraint relations and can build a sparse constraint
    matrix as a :class:`scipy.sparse.coo_matrix`. It also manages groups of
    constrained DOFs, their prescribed values and scaling factors per
    load-case or name.
    """

    def __init__(self, nDofs: int, name: str = "Main") -> None:
        """Initialize a new Constrainer.

        Args:
            nDofs: Total number of DOFs in the global system.
            name: Optional name for the constrainer instance.
        """

        self.nDofs = nDofs
        self.constrainData: Dict[int, List[Any]] = {}
        self.name = name

        # Maps: label -> list of DOF ids
        self.constrainedDofs: Dict[str, List[int]] = {}
        # Maps: label -> list of prescribed values
        self.constrainedVals: Dict[str, List[float]] = {}
        # Maps: label -> scaling factor (per load case)
        self.constrainedFac: Dict[str, float] = {}
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
           
    def addConstraint(self, dofID: int, val: Any, label: str ) -> None:
        """Register a constraint on a DOF.

        The ``val`` argument may be a scalar or a triplet [value, masterID, factor]
        representing a tied DOF. The method updates internal structures that
        are later used to build the constraint matrix.
        """

        if dofID in self.constrainData:
            self.constrainData[dofID].append(val)
        else:
            self.constrainData[dofID] = [val]

        if (type(val) is list) and (len(val) == 3):
            addVal = val[0]
        else:
            addVal = val

        if dofID in self.constrainedDofs[label]:
            self.setFactorForDof(addVal, dofID, label)
            return

        self.constrainedDofs[label].append(dofID)

        self.constrainedVals[label].append(addVal)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def checkConstraints(self, dofspace: Any, nodeTables: Any) -> None:
        """Resolve chained tying relations and expand slave prescriptions.

        This method processes entries in ``self.constrainData`` that indicate
        tied DOFs using the [value, master, factor] convention and flattens
        recursive ties so that slaves with masters that are themselves slaves
        are resolved to base values.
        """

        for item in self.constrainData:

            dofInd = item

            for ilabel in self.constrainedDofs:
                if dofInd in self.constrainedDofs[ilabel]:
                    label = ilabel

            removeItem = []

            for tiedItem in self.constrainData[item]:

                if (type(tiedItem) is list) and len(tiedItem) == 3:
                    Dat = tiedItem

                    valSlave = Dat[0]
                    masterDofID = Dat[1][0]
                    facSlave = Dat[2]

                    if masterDofID in self.constrainData:
                        tempVal: List[float] = []
                        tempFac: List[float] = []

                        # Recursive loop until masterDofID not a list, but prescribed value
                        while masterDofID in self.constrainData:
                            master = self.constrainData[masterDofID][0]
                            if type(master) is list and len(master) == 3:
                                masterDofID = master[1][0]
                                tempVal.append(master[0])
                                tempFac.append(master[2])
                            else:
                                masterFin = master
                                masterDofID = -1

                        for iVal, iFac in reversed(list(zip(tempVal, tempFac))):
                            masterFin += iVal + master * iFac

                        self.addConstraint(dofInd, valSlave + masterFin * facSlave, label)

                        removeItem.append(tiedItem)

            for iRemove in removeItem:
                self.constrainData[item].remove(iRemove)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def flush(self) -> None:
        """Build and store the sparse constraint matrix ``self.C``.

        The returned matrix maps the full DOF vector to the reduced set of
        free DOFs. Slave entries are added to the matrix based on
        previously recorded master/slave relations.
        """

        row: List[int] = []
        col: List[int] = []
        val: List[float] = []
        master: Dict[int, Any] = {}

        iCon = 0

        for iDof in range(self.nDofs):
            if iDof in self.constrainData:
                for item in self.constrainData[iDof]:
                    if type(item) is not list:
                        continue
                    else:
                        if item[1] in self.constrainData:
                            # Something not checked correct in checkConstraint2
                            raise RuntimeError("ERROR - Master of slave is a slave itself")
                        else:
                            master[iDof] = item

            else:
                row.append(iDof)
                col.append(iCon)
                val.append(1.0)

                iCon += 1

        # Assign correct slaves to masters
        for iSlave in master:
            fac = master[iSlave][2]
            masterDofID = master[iSlave][1]

            # Find column for free DOF of the Master
            rowID = row.index(masterDofID)
            col.append(rowID)
            row.append(iSlave)
            val.append(fac)

        self.C = coo_matrix((val, (row, col)), shape=(self.nDofs, iCon))

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def addConstrainedValues(self, a: Any) -> None:
        """Add the constrained values into array ``a`` (incremental).

        The method loops over named constraint groups and increments the
        entries in ``a`` by scaled prescribed values.
        """

        for name in self.constrainedDofs.keys():
            a[self.constrainedDofs[name]] += self.constrainedFac[name] * array(self.constrainedVals[name])

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def setConstrainedValues(self, a: Any) -> None:
        """Overwrite entries in ``a`` with constrained values."""

        for name in self.constrainedDofs.keys():
            a[self.constrainedDofs[name]] = self.constrainedFac[name] * array(self.constrainedVals[name])
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def setConstrainFactor(self, fac: float, loadCase: str = "All_") -> None:
        """Set the constraint scaling factor for all or a specific load case."""

        if loadCase == "All_":
            for name in self.constrainedFac.keys():
                self.constrainedFac[name] = fac
        else:
            self.constrainedFac[loadCase] = fac
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def setPrescribedDofs(self, a: Any, val: float = 0.0) -> None:
        """Set prescribed DOFs in array ``a`` to a scalar value."""

        for name in self.constrainedDofs.keys():
            a[self.constrainedDofs[name]] = array(val)
      
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------      

    def setFactorForDof(self, fac: float, dofID: int, label: str) -> None:
        """Add a factor to the prescribed value for a specific DOF."""

        idx = self.constrainedDofs[label].index(dofID)
        self.constrainedVals[label][idx] += fac
    
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def slaveCount(self) -> int:
        """Return the total number of constrained (slave) DOFs."""

        counter = 0
        for name in self.constrainedFac.keys():
            counter += len(self.constrainedDofs[name])

        return counter