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
from numpy import zeros, ones, ix_, append, repeat, array
from scipy.sparse import coo_matrix
from typing import Any, Tuple, Iterable
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
        rank: Assembly rank (1 for vectors, 2 for matrices).
        action: Name of the element method to call for contributions.

    Returns:
        A tuple containing assembled objects. Types depend on ``rank``.
    """

    # Initialize the global array B with rank 1
    B = zeros(len(globdat.dofs) * ones(1, dtype=int))
    cc = 0.0

    val = array([], dtype=float)
    row = array([], dtype=int)
    col = array([], dtype=int)

    nDof = len(globdat.dofs)

    if action != "commit":
        globdat.resetNodalOutput()

    for elementGroup in globdat.elements.iterGroupNames():
        el_props = getattr(props, elementGroup)

        for iElm, element in enumerate(globdat.elements.iterElementGroup(elementGroup)):
            elemdat = getElementData(element, el_props, globdat)

            elemdat.iElm = iElm
            element.iElm = iElm

            if hasattr(element, "mat"):
                element.mat.reset()

            # Get the element contribution by calling the specified action
            if hasattr(element, action):
                getattr(element, action)(elemdat)

            # Assemble in the global array
            if rank == 1:
                B[elemdat.el_dofs] += elemdat.fint
                cc += elemdat.diss
            elif rank == 2 and action == "getTangentStiffness":
                row = append(row, repeat(elemdat.el_dofs, len(elemdat.el_dofs)))

                for i in range(len(elemdat.el_dofs)):
                    col = append(col, elemdat.el_dofs)

                val = append(val, elemdat.stiff.reshape(len(elemdat.el_dofs) * len(elemdat.el_dofs)))

                B[elemdat.el_dofs] += elemdat.fint
            elif rank == 2 and action == "getMassMatrix":
                row = append(row, repeat(elemdat.el_dofs, len(elemdat.el_dofs)))

                for i in range(len(elemdat.el_dofs)):
                    col = append(col, elemdat.el_dofs)

                val = append(val, elemdat.mass.reshape(len(elemdat.el_dofs) * len(elemdat.el_dofs)))

                B[elemdat.el_dofs] += elemdat.lumped

    globdat.models.run(props, globdat)

    if rank == 1:
        return B, cc
    elif rank == 2:
        return coo_matrix((val, (row, col)), shape=(nDof, nDof)), B


#-------------------------------------------------------------------------------
#  Assemble Internal force
#-------------------------------------------------------------------------------


def assembleInternalForce(props: Properties, globdat: Any) -> Any:
  """Assemble and return the global internal force vector.

  Args:
    props: Global properties container.
    globdat: Global data/state object.

  Returns:
    The assembled internal force array.
  """

  fint = assembleArray(props, globdat, rank=1, action="getInternalForce")
  return fint[0]


#-------------------------------------------------------------------------------
#  Assemble Internal force
#-------------------------------------------------------------------------------


def assembleExternalForce(props: Properties, globdat: Any) -> Any:
    """Assemble and return the global external force vector.

    The external force returned includes contributions assembled from
    elements plus any scaled forcing term stored on ``globdat``.
    """

    fext = assembleArray(props, globdat, rank=1, action="getExternalForce")

    return fext[0] + globdat.fhat * globdat.solverStatus.lam


#-------------------------------------------------------------------------------
#  Assemble Dissipation
#-------------------------------------------------------------------------------
  
  
def assembleDissipation(props: Properties, globdat: Any) -> Any:
    """Assemble and return dissipation contributions.

    Returns the tuple produced by :func:`assembleArray` for the
    ``getDissipation`` action. The first element typically contains
    the assembled vector and the second the scalar dissipation.
    """

    return assembleArray(props, globdat, rank=1, action="getDissipation")
 
 
#-------------------------------------------------------------------------------
#  Assemble Tangent stiffness
#-------------------------------------------------------------------------------


def assembleTangentStiffness(props: Properties, globdat: Any) -> Tuple[Any, Any]:
    """Assemble and return the global tangent stiffness matrix and residual.

    Returns a tuple (stiff_matrix, residual_vector).
    """

    return assembleArray(props, globdat, rank=2, action="getTangentStiffness")


#-------------------------------------------------------------------------------
#  Assemble Mass Matrix
#-------------------------------------------------------------------------------


def assembleMassMatrix(props: Properties, globdat: Any) -> Tuple[Any, Any]:
    """Assemble and return the global mass matrix and lumped vector.

    Returns a tuple (mass_matrix, lumped_vector).
    """

    return assembleArray(props, globdat, rank=2, action="getMassMatrix")


#-------------------------------------------------------------------------------
#  Commit
#-------------------------------------------------------------------------------


def commit(props: Properties, globdat: Any) -> Any:
    """Commit element states by calling the element ``commit`` action.

    Delegates to :func:`assembleArray` with the ``commit`` action.
    """

    return assembleArray(props, globdat, rank=0, action="commit")


#-------------------------------------------------------------------------------
#  getAllConstraints
#-------------------------------------------------------------------------------


def getAllConstraints(props: Properties, globdat: Any) -> None:
    """Invoke `getConstraints` on all elements to collect constraint data.

    Note: this function relies on an ``elemdat`` structure being available
    in the local scope of the caller. The current implementation mirrors the
    original behaviour and intentionally does not change the logic.
    """

    # Loop over the element groups
    for elementGroup in globdat.elements.iterGroupNames():

        # Get the properties corresponding to the elementGroup
        el_props = getattr(props, elementGroup)

        # Pre-create a placeholder element data container so that
        # calls to element.getConstraints can always receive an object.
        elemdat = elementData(np.array([]), np.array([]))

        # Loop over the elements in the elementGroup
        for element in globdat.elements.iterElementGroup(elementGroup):

            # Get the element nodes
            el_nodes = element.getNodes()

            elemdat.nodes = el_nodes
            elemdat.props = el_props

            # Get the element contribution by calling the specified action
            getattr(element, "getConstraints", None)(elemdat)

#-------------------------------------------------------------------------------
#  getElementData
#-------------------------------------------------------------------------------


def getElementData(element: Any, el_props: Properties, globdat: Any) -> elementData:
    """Create and populate an ``elementData`` instance for an element.

    This helper collects node indices, coordinates, DOF indices and the
    current element state from ``globdat`` and stores them in an
    ``elementData`` instance which is returned.
    """

    el_nodes = element.getNodes()

    el_coords = globdat.nodes.getNodeCoords(el_nodes)

    el_dofs = globdat.dofs.getForTypes(el_nodes, element.dofTypes)

    el_a = globdat.state[el_dofs]
    el_Da = globdat.Dstate[el_dofs]

    elemdat = elementData(el_a, el_Da)

    elemdat.coords = el_coords
    elemdat.nodes = el_nodes
    elemdat.props = el_props
    # elemdat.iElm = iElm

    element.globdat = globdat
    # element.iElm = iElm
    elemdat.el_dofs = el_dofs

    if hasattr(element, "matProps"):
        elemdat.matprops = element.matProps

    return elemdat
