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

"""
Bezier shape functions for PyFEM.

This module implements Bezier-based shape functions for finite element analysis,
providing higher-order geometric representations beyond standard polynomial shape
functions. Bezier elements are useful for modeling curved geometries and achieving
geometric exactness in isogeometric analysis.

Module contents:
    - getBezierLine4: Compute Bezier shape functions for 4-node cubic Bezier curve
    - calcWeight: Compute integration weight from Jacobian
    - getElemBezierData: Generate shape data for Bezier elements with integration

References:
    Isogeometric analysis and related computational mechanics applications.
"""

from typing import Union, Tuple
import numpy as np
from scipy.linalg import norm, det, inv

try:
    # For newer versions of scipy
    from scipy.special import p_roots as gauss_scheme
except ImportError:
    # For older versions of scipy
    from scipy.special.orthogonal import p_roots as gauss_scheme

from .shapeFunctions import shapeData, elemShapeData, getIntegrationPoints, getElemType


def getBezierLine4(xi: Union[float, np.floating], C: np.ndarray) -> shapeData:
    """
    Compute shape functions and derivatives for a cubic Bezier line element.

    Evaluates the shape functions and their derivatives with respect to the
    parametric coordinate xi for a 4-node cubic Bezier curve. The physical
    coordinates are computed using the control points in matrix C.

    Args:
        xi (float): Parametric coordinate in the range [-1, 1].
        C (np.ndarray): Control points matrix of shape (nDim, 4), where each column
            is a control point in physical space. nDim is typically 2 or 3.

    Returns:
        shapeData: Object containing:
            - h: Physical coordinates at parameter xi (shape: (nDim,))
            - dhdxi: Derivative of physical coordinates w.r.t. xi (shape: (nDim, 1))
            - xi: The input parametric coordinate

    Raises:
        NotImplementedError: If xi is a 1D array or C does not have exactly 4 columns.

    Notes:
        The cubic Bezier basis functions are:
        - B[0] = -0.125*(xi-1)^3
        - B[1] = 0.375*(xi-1)^2*(xi+1)
        - B[2] = -0.375*(xi-1)*(xi+1)^2
        - B[3] = 0.125*(xi+1)^3

    Examples:
        >>> import numpy as np
        >>> xi = 0.0  # Evaluate at parametric center
        >>> C = np.array([[0, 1, 2, 3], [0, 1, 1, 0]])  # 2D control points
        >>> sData = getBezierLine4(xi, C)
        >>> sData.h.shape
        (2,)
    """
    # Check the dimensions of the parametric space
    if type(xi) == 'numpy.float64':
        raise NotImplementedError('1D only')
    if C.shape[1] != 4:
        raise NotImplementedError('C needs to have 4 columns.')

    sData = shapeData()

    # Store parametric coordinate
    sData.xi = xi

    # Initialize basis function and derivative arrays
    B = np.empty(4)
    dBdxi = np.empty(shape=(4, 1))

    # Calculate cubic Bezier basis functions
    B[0] = -0.125 * (xi - 1.) ** 3
    B[1] = 0.375 * (xi - 1.) ** 2 * (xi + 1.)
    B[2] = -0.375 * (xi - 1.) * (xi + 1.) ** 2
    B[3] = 0.125 * (xi + 1.) ** 3

    # Calculate derivatives of basis functions w.r.t. parametric coordinate
    dBdxi[0, 0] = -0.375 * (xi - 1.) ** 2
    dBdxi[1, 0] = 0.75 * (xi - 1.0) * (xi + 1.0) + 0.375 * (xi - 1.) ** 2
    dBdxi[2, 0] = -0.375 * (1 + xi) ** 2 - 0.75 * (1 + xi) * (xi - 1)
    dBdxi[3, 0] = 0.375 * (xi + 1.) ** 2

    # Map from parametric space to physical space
    sData.h = C @ B
    sData.dhdxi = C @ dBdxi

    return sData


def calcWeight(jac: np.ndarray) -> Union[np.floating, float]:
    """
    Calculate integration weight from Jacobian matrix.

    Computes the determinant of the Jacobian for volume integration (square matrix)
    or the norm for line integration (rectangular matrix), which are used to weight
    integration points in numerical quadrature.

    Args:
        jac (np.ndarray): Jacobian matrix. Can be:
            - Square matrix (nDim, nDim) for volume integration: returns determinant
            - Rectangular matrix (1, 2) for 1D curve in 2D space: returns curve length

    Returns:
        np.floating | float: The integration weight. For square Jacobians, returns
            the determinant. For rectangular matrices, returns the norm (magnitude)
            of the Jacobian.

    Notes:
        The weight is typically multiplied by integration point weights from
        quadrature rules to compute integrals over curved elements.

    Examples:
        >>> import numpy as np
        >>> # 2D Jacobian (square)
        >>> jac_2d = np.array([[1.0, 0.0], [0.0, 1.0]])
        >>> w = calcWeight(jac_2d)
        >>> w
        1.0

        >>> # 1D curve in 2D space
        >>> jac_1d = np.array([[1.0, 1.0]])
        >>> w = calcWeight(jac_1d)
        >>> np.isclose(w, np.sqrt(2.0))
        True
    """
    n = jac.shape

    if n[0] == n[1]:
        # Square Jacobian: return determinant for volume integration
        return det(jac)
    elif n[0] == 1 and n[1] == 2:
        # 1D curve in 2D: return norm (arc length differential)
        return np.sqrt(np.sum(np.sum(jac * jac)))


def getElemBezierData(
    elemCoords: np.ndarray,
    C: np.ndarray,
    order: int = 4,
    method: str = "Gauss",
    elemType: str = 'default'
) -> elemShapeData:
    """
    Generate shape data and integration information for Bezier elements.

    Evaluates shape functions and their spatial derivatives at integration points
    for a Bezier element. This function is used to assemble element stiffness
    matrices and load vectors in isogeometric analysis.

    Args:
        elemCoords (np.ndarray): Physical coordinates of element nodes, shape
            (nDim, nNodes). These are the physical nodes in the background mesh.
        C (np.ndarray): Control points matrix, shape (nDim, nControlPts).
            Defines the geometric mapping via Bezier basis functions.
        order (int, optional): Integration order (number of integration points).
            Defaults to 4. Higher orders give more accurate integration.
        method (str, optional): Quadrature rule. Defaults to "Gauss" for
            Gauss-Legendre quadrature.
        elemType (str, optional): Element type identifier (e.g., 'Line4', 'Line3').
            If 'default', automatically detected from elemCoords shape.
            Defaults to 'default'.

    Returns:
        elemShapeData: Container with integration point data. Attributes:
            - sData: List of shapeData objects, one per integration point, each
              containing:
              - h: Physical coordinates at integration point (shape: (nDim,))
              - dhdxi: Derivatives w.r.t. parametric coord (shape: (nDim, 1))
              - dhdx: Derivatives w.r.t. physical coords (shape: (nDim, nDim))
              - weight: Integration weight including Jacobian determinant
              - xi: Parametric coordinate of integration point

    Raises:
        NotImplementedError: If elemType is unknown or not supported.

    Notes:
        Integration points are obtained from Line3 element quadrature rules.
        The Jacobian is computed from the physical derivatives of shape functions.
        The spatial derivatives (dhdx) are computed only for square Jacobians.

    Examples:
        >>> import numpy as np
        >>> # Control points for a curved quadratic line
        >>> C = np.array([[0, 0.5, 1.0], [0, 0.1, 0.0]])
        >>> elemCoords = np.array([[0, 1.0], [0, 0.0]])
        >>> elemData = getElemBezierData(elemCoords, C, order=2, elemType='Line2')
        >>> len(elemData.sData)  # Number of integration points
        2
        >>> elemData.sData[0].weight > 0  # Positive weight
        True
    """
    elemData = elemShapeData()

    if elemType == 'default':
        elemType = getElemType(elemCoords)

    (intCrds, intWghts) = getIntegrationPoints("Line3", order, method)

    for xi, intWeight in zip(np.real(intCrds), intWghts):
        try:
            sData = eval('getBezier' + elemType + '(xi,C)')
        except:
            raise NotImplementedError('Unknown type :' + elemType)

        jac = sData.dhdxi.T @ elemCoords

        if jac.shape[0] is jac.shape[1]:
            sData.dhdx = (inv(jac) @ sData.dhdxi.T).T

        sData.weight = calcWeight(jac) * intWeight

        elemData.sData.append(sData)

    return elemData

