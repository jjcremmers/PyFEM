# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

import numpy as np
from numpy import zeros, ones, ix_, append, repeat, array
from scipy.sparse import coo_matrix
from typing import Any, Tuple, Iterable, Optional
from numpy.typing import NDArray
from pyfem.util.dataStructures import Properties
from pyfem.util.dataStructures import elementData


#-------------------------------------------------------------------------------
#  Matrix builder
#-------------------------------------------------------------------------------


class ModelBuilder:
    """
    Incrementally build sparse matrix data in COO format for finite element assembly.
    Also stores a global vector B and a scalar c for additional model data.

    Attributes:
        nDofs (int): Number of global degrees of freedom.
        val (np.ndarray): Values of the matrix entries.
        row (np.ndarray): Row indices for COO format.
        col (np.ndarray): Column indices for COO format.
        B (np.ndarray): Global vector (initialized to zeros).
        c (float): Scalar value (initialized to 0.0).
    """

    def __init__(self, nDofs: int) -> None:
        """
        Initialize a new builder for a system with nDofs.

        Args:
            nDofs (int): Number of global degrees of freedom.
        """
        self.nDofs: int = nDofs
        self.clear()

    def clear(self) -> None:
        """
        Reset stored row/col/value arrays, B, and c.
        """
        self.val: NDArray[np.floating] = array([], dtype=float)
        self.row: NDArray[np.integer] = array([], dtype=int)
        self.col: NDArray[np.integer] = array([], dtype=int)
        self.B: NDArray[np.floating] = zeros(self.nDofs, dtype=float)
        self.c: float = 0.0

    def append(self, a: NDArray[np.floating], dofs: NDArray[np.integer]) -> None:
        """
        Append element matrix `a` using associated dof indices.

        Args:
            a (np.ndarray): Element matrix (square, shape [n, n]).
            dofs (np.ndarray): DOF indices for the element (length n).
        """
        n = len(dofs)
        self.row = append(self.row, repeat(dofs, n))
        self.col = append(self.col, np.tile(dofs, n))
        self.val = append(self.val, a.reshape(n * n))

    def getMatrix(self) -> coo_matrix:
        """
        Return the assembled COO sparse matrix.

        Returns:
            scipy.sparse.coo_matrix: Assembled sparse matrix.
        """
        return coo_matrix((self.val, (self.row, self.col)), shape=(self.nDofs, self.nDofs))


#-------------------------------------------------------------------------------
#  Assemble Internal force
#-------------------------------------------------------------------------------


def assembleInternalForce(props: Properties, globdat: Any) -> NDArray[np.floating]:
    """
    Assemble and return the global internal force vector.

    Computes the internal force vector by calling the 'getInternalForce'
    method on all elements and assembling their contributions.

    Args:
        props (Properties): Global properties container.
        globdat (Any): Global data/state object.

    Returns:
        np.ndarray: The assembled internal force vector.
    """

    mbuilder = ModelBuilder(len(globdat.dofs))

    globdat.resetNodalOutput()

    for elementGroup in globdat.elements.iterGroupNames():
        el_props = getattr(props, elementGroup)

        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            elemdat = getElementData(iElm, element, el_props, globdat)

            if hasattr(element, "mat"):
                element.mat.reset()

            if hasattr(element, "getInternalForce"):
                element.getInternalForce(elemdat)

            mbuilder.B[elemdat.el_dofs] += elemdat.fint

    globdat.models.takeAction( "getInternalForce" , mbuilder , props , globdat )

    return mbuilder.B


#-------------------------------------------------------------------------------
#  Assemble Internal force
#-------------------------------------------------------------------------------


def assembleExternalForce(props: Properties, globdat: Any) -> NDArray[np.floating]:
    """
    Assemble and return the global external force vector.

    Computes the external force vector by calling the 'getExternalForce'
    method on all elements and assembling their contributions. The external
    force returned includes contributions assembled from elements plus any
    scaled forcing term stored on ``globdat``.

    Args:
        props (Properties): Global properties container.
        globdat (Any): Global data/state object.

    Returns:
        np.ndarray: The assembled external force vector, including
        the scaled load factor contribution (globdat.fhat * globdat.solverStatus.lam).
    """

    mbuilder = ModelBuilder(len(globdat.dofs))

    globdat.resetNodalOutput()

    for elementGroup in globdat.elements.iterGroupNames():
        el_props = getattr(props, elementGroup)

        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            elemdat = getElementData(iElm, element, el_props, globdat)

            if hasattr(element, "mat"):
                element.mat.reset()

            if hasattr(element, "getExternalForce"):
                element.getExternalForce(elemdat)

            mbuilder.B[elemdat.el_dofs] += elemdat.fint

    globdat.models.takeAction( "getExternalForce" , mbuilder , props , globdat )

    return mbuilder.B + globdat.fhat * globdat.solverStatus.lam


#-------------------------------------------------------------------------------
#  Assemble Dissipation
#-------------------------------------------------------------------------------
  
  
def assembleDissipation(props: Properties, globdat: Any) -> Tuple[NDArray[np.floating], float]:
    """
    Assemble and return dissipation contributions.

    Computes dissipation by calling the 'getDissipation' method on all elements.

    Args:
        props (Properties): Global properties container.
        globdat (Any): Global data/state object.

    Returns:
        tuple[np.ndarray, float]:
            - dissipation_vector: Assembled dissipation force vector
            - accumulated_dissipation: Total scalar dissipation from all elements
    """

    mbuilder = ModelBuilder(len(globdat.dofs))

    globdat.resetNodalOutput()

    for elementGroup in globdat.elements.iterGroupNames():
        el_props = getattr(props, elementGroup)

        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            elemdat = getElementData(iElm, element, el_props, globdat)

            if hasattr(element, "mat"):
                element.mat.reset()

            if hasattr(element, "getDissipation"):
                element.getDissipation(elemdat)

            mbuilder.B[elemdat.el_dofs] += elemdat.fint
            mbuilder.c += elemdat.diss

    globdat.models.takeAction( "getDissipation" , mbuilder , props , globdat )

    return mbuilder.B, mbuilder.c
 
 
#-------------------------------------------------------------------------------
#  Assemble Tangent stiffness
#-------------------------------------------------------------------------------


def assembleTangentStiffness(props: Properties, globdat: Any) -> Tuple[coo_matrix, NDArray[np.floating]]:
    """
    Assemble and return the global tangent stiffness matrix and residual.

    Computes the tangent stiffness matrix by calling the 'getTangentStiffness'
    method on all elements and assembling their contributions into a sparse matrix.

    Args:
        props (Properties): Global properties container.
        globdat (Any): Global data/state object.

    Returns:
        tuple[scipy.sparse.coo_matrix, np.ndarray]:
            - stiff_matrix: Global tangent stiffness matrix in COO sparse format
            - residual_vector: Assembled internal force residual vector
    """

    mbuilder = ModelBuilder(len(globdat.dofs))

    globdat.resetNodalOutput()

    for elementGroup in globdat.elements.iterGroupNames():
        el_props = getattr(props, elementGroup)

        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            elemdat = getElementData(iElm, element, el_props, globdat)
 
            if hasattr(element, "mat"):
                element.mat.reset()

            if hasattr(element, "getTangentStiffness"):
                element.getTangentStiffness(elemdat)

            mbuilder.append(elemdat.stiff, elemdat.el_dofs)

            mbuilder.B[elemdat.el_dofs] += elemdat.fint
    
    globdat.models.takeAction( "getTangentStiffness" , mbuilder , props , globdat )

    return mbuilder.getMatrix(), mbuilder.B


#-------------------------------------------------------------------------------
#  Assemble Mass Matrix
#-------------------------------------------------------------------------------


def assembleMassMatrix(props: Properties, globdat: Any) -> Tuple[coo_matrix, NDArray[np.floating]]:
    """
    Assemble and return the global mass matrix and lumped mass vector.

    Computes the mass matrix by calling the 'getMassMatrix' method on all
    elements and assembling their contributions into a sparse matrix.

    Args:
        props (Properties): Global properties container.
        globdat (Any): Global data/state object.

    Returns:
        tuple[scipy.sparse.coo_matrix, np.ndarray]:
            - mass_matrix: Global mass matrix in COO sparse format
            - lumped_mass_vector: Assembled lumped mass vector (diagonal approximation)
    """

    mbuilder = ModelBuilder(len(globdat.dofs))

    globdat.resetNodalOutput()

    for elementGroup in globdat.elements.iterGroupNames():
        el_props = getattr(props, elementGroup)

        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            elemdat = getElementData(iElm, element, el_props, globdat)

            if hasattr(element, "mat"):
                element.mat.reset()

            if hasattr(element, "getMassMatrix"):
                element.getMassMatrix(elemdat)

            mbuilder.append(elemdat.mass, elemdat.el_dofs)

            mbuilder.B[elemdat.el_dofs] += elemdat.lumped
    
    globdat.models.takeAction( "getMassMatrix" , mbuilder , props , globdat )

    return mbuilder.getMatrix(), mbuilder.B

#-------------------------------------------------------------------------------
#  Commit
#-------------------------------------------------------------------------------


def commit(props: Properties, globdat: Any) -> None:
    """
    Commit element states by calling the element 'commit' method.

    This function is called after a successful time step or load step to
    finalize and store the current element states (e.g., history variables,
    plastic strains, damage parameters).

    Args:
        props (Properties): Global properties container.
        globdat (Any): Global data/state object.
    """

    mbuilder = ModelBuilder(len(globdat.dofs))
    
    for elementGroup in globdat.elements.iterGroupNames():
        el_props = getattr(props, elementGroup)

        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            elemdat = getElementData(iElm, element, el_props, globdat)

            if hasattr(element, "mat"):
                element.mat.reset()

            if hasattr(element, "commit"):
                element.commit(elemdat)

    globdat.models.takeAction( "commit" , mbuilder , props , globdat )

    return None


#-------------------------------------------------------------------------------
#  getAllConstraints
#-------------------------------------------------------------------------------


def getAllConstraints(props: Properties, globdat: Any) -> None:
    """
    Invoke 'getConstraints' on all elements to collect constraint data.

    This function iterates over all element groups and elements, calling their
    'getConstraints' method if available. This is typically used for multi-point
    constraints, contact constraints, or other element-level constraint definitions.

    Args:
        props (Properties): Global properties container.
        globdat (Any): Global data/state object.

    Note:
        The current implementation creates a minimal elemdat structure for each
        element. This mirrors the original behavior and intentionally does not
        change the logic.
    """

    # Loop over all element groups
    for elementGroup in globdat.elements.iterGroupNames():

        # Get the properties corresponding to the element group
        el_props = getattr(props, elementGroup)

        # Pre-create a placeholder element data container so that
        # calls to element.getConstraints can always receive an object
        elemdat = elementData(np.array([]), np.array([]))

        # Loop over all elements in the element group
        for element in globdat.elements.iterElementGroup(elementGroup):

            # Get the element node indices
            el_nodes = element.getNodes()

            # Populate element data with nodes and properties
            elemdat.nodes = el_nodes
            elemdat.props = el_props

            # Call the getConstraints method if it exists on the element
            getattr(element, "getConstraints", None)(elemdat)

#-------------------------------------------------------------------------------
#  getElementData
#-------------------------------------------------------------------------------


def getElementData(iElm: int, element: Any, el_props: Properties, globdat: Any) -> elementData:
    """
    Create and populate an elementData instance for an element.

    This helper function gathers all necessary data for an element from the
    global data structure, including:
        - Node indices and coordinates
        - Degree of freedom indices and values
        - Current state vector and state increment
        - Element properties and material properties

    Args:
        iElm (int): Element index in the group.
        element (Any): The element object for which to gather data.
        el_props (Properties): Properties object for the element's group.
        globdat (Any): Global data/state object containing mesh, DOFs, and state.

    Returns:
        elementData: An instance populated with all element-specific information
        needed for element computations.
    """

    # Get element node indices
    el_nodes = element.getNodes()

    # Get element node coordinates from global node array
    el_coords = globdat.nodes.getNodeCoords(el_nodes)

    # Get element DOF indices based on node indices and DOF types
    el_dofs = globdat.dofs.getForTypes(el_nodes, element.dofTypes)

    # Extract current state and state increment for element DOFs
    el_a = globdat.state[el_dofs]
    el_Da = globdat.Dstate[el_dofs]

    # Create elementData object with state vectors
    elemdat = elementData(el_a, el_Da)

    # Populate elementData with additional information
    elemdat.coords = el_coords
    elemdat.nodes = el_nodes
    elemdat.props = el_props
    elemdat.el_dofs = el_dofs
    elemdat.iElm    = iElm

    # Attach global data to element for access if needed
    element.globdat = globdat

    # Add material properties if element has them
    if hasattr(element, "matProps"):
        elemdat.matprops = element.matProps

    return elemdat
