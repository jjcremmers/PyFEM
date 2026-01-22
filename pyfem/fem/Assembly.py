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
#  Assemble Internal force
#-------------------------------------------------------------------------------

def assembleArray(props: Properties, globdat: Any, rank: int, action: str) -> Tuple[Any, Any]:
    """Assemble global arrays (forces or matrices) from element contributions.

    This routine iterates over element groups and elements, calls the
    requested element action (e.g. ``getInternalForce``, ``getTangentStiffness``)
    and assembles the returned element-level vectors or matrices into
    global arrays. The function returns a tuple; for vector assembly the
    tuple is (B, cc) where B is the assembled array and cc the accumulated
    dissipation. For matrix assembly the tuple is (sparse_matrix, B) where
    B is the assembled residual vector.

    Args:
        props: Global properties container (typically a ``Properties`` instance).
        globdat: Global data/state object used by elements.
        rank: Assembly rank (1 for vectors, 2 for matrices, 0 for commit).
        action: Name of the element method to call for contributions.
               Common actions include:
               - 'getInternalForce': Compute internal force contributions
               - 'getExternalForce': Compute external force contributions
               - 'getTangentStiffness': Compute tangent stiffness matrix
               - 'getMassMatrix': Compute mass matrix
               - 'getDissipation': Compute dissipation contributions
               - 'commit': Commit element states

    Returns:
        A tuple containing assembled objects. Types depend on ``rank``:
        - rank == 1: (assembled_vector, accumulated_dissipation)
        - rank == 2: (sparse_matrix, residual_vector)
        - rank == 0: (None, None) for commit operations
    """

    # Initialize the global array B (assembled vector)
    B = zeros(len(globdat.dofs) * ones(1, dtype=int))
    
    # Initialize accumulated dissipation
    cc = 0.0

    # Initialize sparse matrix storage arrays (COO format)
    val = array([], dtype=float)
    row = array([], dtype=int)
    col = array([], dtype=int)

    # Store total number of degrees of freedom
    nDof = len(globdat.dofs)

    # Reset nodal output for all actions except commit
    if action != "commit":
        globdat.resetNodalOutput()

    # Loop over all element groups
    for elementGroup in globdat.elements.iterGroupNames():
        # Get properties for current element group
        el_props = getattr(props, elementGroup)

        # Loop over all elements in the current group
        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            # Get element data (nodes, DOFs, coordinates, state)
            elemdat = getElementData(element, el_props, globdat)

            # Store element index
            elemdat.iElm = iElm
            element.iElm = iElm

            # Reset material state if element has a material
            if hasattr(element, "mat"):
                element.mat.reset()

            # Get the element contribution by calling the specified action
            if hasattr(element, action):
                getattr(element, action)(elemdat)

            # Assemble element contributions into global arrays
            if rank == 1:
                # Vector assembly: add internal forces and accumulate dissipation
                B[elemdat.el_dofs] += elemdat.fint
                cc += elemdat.diss
            elif rank == 2 and action == "getTangentStiffness":
                # Matrix assembly: store stiffness matrix in COO format
                row = append(row, repeat(elemdat.el_dofs, len(elemdat.el_dofs)))

                for i in range(len(elemdat.el_dofs)):
                    col = append(col, elemdat.el_dofs)

                val = append(val, elemdat.stiff.reshape(len(elemdat.el_dofs) * len(elemdat.el_dofs)))

                # Also assemble internal force residual
                B[elemdat.el_dofs] += elemdat.fint
            elif rank == 2 and action == "getMassMatrix":
                # Matrix assembly: store mass matrix in COO format
                row = append(row, repeat(elemdat.el_dofs, len(elemdat.el_dofs)))

                for i in range(len(elemdat.el_dofs)):
                    col = append(col, elemdat.el_dofs)

                val = append(val, elemdat.mass.reshape(len(elemdat.el_dofs) * len(elemdat.el_dofs)))

                # Also assemble lumped mass vector
                B[elemdat.el_dofs] += elemdat.lumped

    # Run any additional models (constraints, boundary conditions, etc.)

    globdat.B = B
    globdat.val = val
    globdat.row = row
    globdat.col = col

    print("before ", B)

    globdat.models.run(props, globdat)

    print("after ", B)


    # Return appropriate result based on assembly rank
    if rank == 1:
        # Vector assembly: return assembled vector and dissipation
        return B, cc
    elif rank == 2:
        # Matrix assembly: return sparse matrix and residual vector
        return coo_matrix((val, (row, col)), shape=(nDof, nDof)), B


#-------------------------------------------------------------------------------
#  Assemble Internal force
#-------------------------------------------------------------------------------


def assembleInternalForce(props: Properties, globdat: Any) -> NDArray[np.floating]:
    """Assemble and return the global internal force vector.

    Computes the internal force vector by calling the 'getInternalForce'
    action on all elements and assembling their contributions.

    Args:
        props: Global properties container.
        globdat: Global data/state object.

    Returns:
        The assembled internal force vector as a numpy array.
    """

    fint = assembleArray(props, globdat, rank=1, action="getInternalForce")
    return fint[0]


#-------------------------------------------------------------------------------
#  Assemble Internal force
#-------------------------------------------------------------------------------


def assembleExternalForce(props: Properties, globdat: Any) -> NDArray[np.floating]:
    """Assemble and return the global external force vector.

    Computes the external force vector by calling the 'getExternalForce'
    action on all elements and assembling their contributions. The external
    force returned includes contributions assembled from elements plus any
    scaled forcing term stored on ``globdat``.

    Args:
        props: Global properties container.
        globdat: Global data/state object.

    Returns:
        The assembled external force vector as a numpy array, including
        the scaled load factor contribution (globdat.fhat * globdat.solverStatus.lam).
    """

    fext = assembleArray(props, globdat, rank=1, action="getExternalForce")

    return fext[0] + globdat.fhat * globdat.solverStatus.lam


#-------------------------------------------------------------------------------
#  Assemble Dissipation
#-------------------------------------------------------------------------------
  
  
def assembleDissipation(props: Properties, globdat: Any) -> Tuple[NDArray[np.floating], float]:
    """Assemble and return dissipation contributions.

    Computes dissipation by calling the 'getDissipation' action on all elements.

    Args:
        props: Global properties container.
        globdat: Global data/state object.

    Returns:
        A tuple (dissipation_vector, accumulated_dissipation) where:
        - dissipation_vector: Assembled dissipation force vector
        - accumulated_dissipation: Total scalar dissipation from all elements
    """

    return assembleArray(props, globdat, rank=1, action="getDissipation")
 
 
#-------------------------------------------------------------------------------
#  Assemble Tangent stiffness
#-------------------------------------------------------------------------------


def assembleTangentStiffness(props: Properties, globdat: Any) -> Tuple[Any, NDArray[np.floating]]:
    """Assemble and return the global tangent stiffness matrix and residual.

    Computes the tangent stiffness matrix by calling the 'getTangentStiffness'
    action on all elements and assembling their contributions into a sparse matrix.

    Args:
        props: Global properties container.
        globdat: Global data/state object.

    Returns:
        A tuple (stiff_matrix, residual_vector) where:
        - stiff_matrix: Global tangent stiffness matrix in COO sparse format
        - residual_vector: Assembled internal force residual vector
    """

    return assembleArray(props, globdat, rank=2, action="getTangentStiffness")


#-------------------------------------------------------------------------------
#  Assemble Mass Matrix
#-------------------------------------------------------------------------------


def assembleMassMatrix(props: Properties, globdat: Any) -> Tuple[Any, NDArray[np.floating]]:
    """Assemble and return the global mass matrix and lumped mass vector.

    Computes the mass matrix by calling the 'getMassMatrix' action on all
    elements and assembling their contributions into a sparse matrix.

    Args:
        props: Global properties container.
        globdat: Global data/state object.

    Returns:
        A tuple (mass_matrix, lumped_mass_vector) where:
        - mass_matrix: Global mass matrix in COO sparse format
        - lumped_mass_vector: Assembled lumped mass vector (diagonal approximation)
    """

    return assembleArray(props, globdat, rank=2, action="getMassMatrix")


#-------------------------------------------------------------------------------
#  Commit
#-------------------------------------------------------------------------------


def commit(props: Properties, globdat: Any) -> Tuple[Any, Any]:
    """Commit element states by calling the element 'commit' action.

    This function is called after a successful time step or load step to
    finalize and store the current element states (e.g., history variables,
    plastic strains, damage parameters).

    Args:
        props: Global properties container.
        globdat: Global data/state object.

    Returns:
        A tuple returned by assembleArray (typically unused for commit operations).
    """

    return assembleArray(props, globdat, rank=0, action="commit")


#-------------------------------------------------------------------------------
#  getAllConstraints
#-------------------------------------------------------------------------------


def getAllConstraints(props: Properties, globdat: Any) -> None:
    """Invoke 'getConstraints' on all elements to collect constraint data.

    This function iterates over all element groups and elements, calling their
    'getConstraints' method if available. This is typically used for multi-point
    constraints, contact constraints, or other element-level constraint definitions.

    Args:
        props: Global properties container.
        globdat: Global data/state object.

    Returns:
        None. Constraints are typically added to globdat during element calls.

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


def getElementData(element: Any, el_props: Properties, globdat: Any) -> elementData:
    """Create and populate an elementData instance for an element.

    This helper function gathers all necessary data for an element from the
    global data structure, including:
    - Node indices and coordinates
    - Degree of freedom indices and values
    - Current state vector and state increment
    - Element properties and material properties

    Args:
        element: The element object for which to gather data.
        el_props: Properties object for the element's group.
        globdat: Global data/state object containing mesh, DOFs, and state.

    Returns:
        An elementData instance populated with all element-specific information
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

    # Attach global data to element for access if needed
    element.globdat = globdat

    # Add material properties if element has them
    if hasattr(element, "matProps"):
        elemdat.matprops = element.matProps

    return elemdat
