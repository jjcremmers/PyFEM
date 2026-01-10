# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from math import sqrt
import numpy as np
from scipy.linalg import norm, det, inv

try:
    # For newer versions of scipy
    from scipy.special import p_roots as gauss_scheme
except ImportError:
    # For older versions of scipy
    from scipy.special.orthogonal import p_roots as gauss_scheme

from typing import List, Tuple


class shapeData:
    """
    Container for shape function data at a single integration point.
    
    Attributes:
        h: Shape function values at the integration point
        dhdxi: Derivatives of shape functions w.r.t. parent coordinates
        dhdx: Derivatives of shape functions w.r.t. physical coordinates
        xi: Coordinates in parent element
        x: Coordinates in physical element
        weight: Integration weight
    """
    pass


class elemShapeData:
    """
    Container for shape function data for an entire element.
    
    Stores shape data for all integration points and provides iteration
    capability over the integration points.
    """

    def __init__(self) -> None:
        """Initialize empty list of shape data."""
        self.sData: List[shapeData] = []

    def __iter__(self):
        """Iterate over integration points in the element."""
        return iter(self.sData)

    def __len__(self) -> int:
        """Return the number of integration points in the element."""
        return len(self.sData)


#----------------------------------------------------------------------

def getShapeLine2(xi: float) -> shapeData:
    """
    Shape functions for 1D linear line element with 2 nodes (Line2).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of the integration point in parent element
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not a float (1D only)
    """

    if not isinstance(xi, (float, np.floating)):
        raise NotImplementedError(f'1D only: xi is of type {type(xi)}')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(2)
    sData.dhdxi = np.empty(shape=(2, 1))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 0.5 * (1.0 - xi)
    sData.h[1] = 0.5 * (1.0 + xi)

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -0.5
    sData.dhdxi[1, 0] = 0.5

    return sData

#----------------------------------------------------------------------

def getShapeLine3(xi: float) -> shapeData:
    """
    Shape functions for 1D quadratic line element with 3 nodes (Line3).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of the integration point in parent element
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not a float (1D only)
    """

    if not isinstance(xi, (float, np.floating)):
        raise NotImplementedError(f'1D only: xi is of type {type(xi)}')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(3)
    sData.dhdxi = np.empty(shape=(1, 3))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 0.5 * (1.0 - xi) - 0.5 * (1.0 - xi * xi)
    sData.h[1] = 1 - xi * xi
    sData.h[2] = 0.5 * (1.0 + xi) - 0.5 * (1.0 - xi * xi)

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -0.5 + xi
    sData.dhdxi[0, 1] = -2.0 * xi
    sData.dhdxi[0, 2] = 0.5 + xi

    return sData

#----------------------------------------------------------------------

def getShapeTria3(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 2D triangular element with 3 nodes (Tria3).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of integration point as ndarray of length 2
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 2D (length != 2)
    """

    if len(xi) != 2:
        raise NotImplementedError('2D only')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(3)
    sData.dhdxi = np.empty(shape=(3, 2))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 1.0 - xi[0] - xi[1]
    sData.h[1] = xi[0]
    sData.h[2] = xi[1]

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -1.0
    sData.dhdxi[1, 0] = 1.0
    sData.dhdxi[2, 0] = 0.0

    sData.dhdxi[0, 1] = -1.0
    sData.dhdxi[1, 1] = 0.0
    sData.dhdxi[2, 1] = 1.0

    return sData

#-------------------------------------

def getShapeQuad4(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 2D quadrilateral element with 4 nodes (Quad4).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of integration point as ndarray of length 2
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 2D (length != 2)
    """

    if len(xi) != 2:
        raise NotImplementedError('2D only')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(4)
    sData.dhdxi = np.empty(shape=(4, 2))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 0.25 * (1.0 - xi[0]) * (1.0 - xi[1])
    sData.h[1] = 0.25 * (1.0 + xi[0]) * (1.0 - xi[1])
    sData.h[2] = 0.25 * (1.0 + xi[0]) * (1.0 + xi[1])
    sData.h[3] = 0.25 * (1.0 - xi[0]) * (1.0 + xi[1])

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -0.25 * (1.0 - xi[1])
    sData.dhdxi[1, 0] = 0.25 * (1.0 - xi[1])
    sData.dhdxi[2, 0] = 0.25 * (1.0 + xi[1])
    sData.dhdxi[3, 0] = -0.25 * (1.0 + xi[1])

    sData.dhdxi[0, 1] = -0.25 * (1.0 - xi[0])
    sData.dhdxi[1, 1] = -0.25 * (1.0 + xi[0])
    sData.dhdxi[2, 1] = 0.25 * (1.0 + xi[0])
    sData.dhdxi[3, 1] = 0.25 * (1.0 - xi[0])

    return sData

#-------------------------------------

def getShapeTria6(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 2D triangular element with 6 nodes (Tria6).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of integration point as ndarray of length 2
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 2D (length != 2)
    """

    if len(xi) != 2:
        raise NotImplementedError('2D only')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(6)
    sData.dhdxi = np.empty(shape=(6, 2))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 1.0 - xi[0] - xi[1] - 2.0 * xi[0] * (1.0 - xi[0] - xi[1]) - 2.0 * xi[1] * (1.0 - xi[0] - xi[1])
    sData.h[1] = xi[0] - 2.0 * xi[0] * (1.0 - xi[0] - xi[1]) - 2.0 * xi[0] * xi[1]
    sData.h[2] = xi[1] - 2.0 * xi[0] * xi[1] - 2.0 * xi[1] * (1.0 - xi[0] - xi[1])
    sData.h[3] = 4.0 * xi[0] * (1.0 - xi[0] - xi[1])
    sData.h[4] = 4.0 * xi[0] * xi[1]
    sData.h[5] = 4.0 * xi[1] * (1.0 - xi[0] - xi[1])

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -1.0 - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[0] + 2.0 * xi[1]
    sData.dhdxi[1, 0] = 1.0 - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[0] - 2.0 * xi[1]
    sData.dhdxi[2, 0] = 0.0
    sData.dhdxi[3, 0] = 4.0 * (1.0 - xi[0] - xi[1]) - 4.0 * xi[0]
    sData.dhdxi[4, 0] = 4.0 * xi[1]
    sData.dhdxi[5, 0] = -4.0 * xi[1]

    sData.dhdxi[0, 1] = -1.0 + 2.0 * xi[0] - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[1]
    sData.dhdxi[1, 1] = 0.0
    sData.dhdxi[2, 1] = 1.0 - 2.0 * xi[0] - 2.0 * (1.0 - xi[0] - xi[1]) + 2.0 * xi[1]
    sData.dhdxi[3, 1] = -4.0 * xi[0]
    sData.dhdxi[4, 1] = 4.0 * xi[0]
    sData.dhdxi[5, 1] = 4.0 * (1.0 - xi[0] - xi[1]) - 4.0 * xi[1]

    return sData

#-------------------------------------

def getShapeQuad8(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 2D quadrilateral element with 8 nodes (Quad8).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of integration point as ndarray of length 2
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 2D (length != 2)
    """

    if len(xi) != 2:
        raise NotImplementedError('2D only')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(8)
    sData.dhdxi = np.empty(shape=(8, 2))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = -0.25 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[0] + xi[1])
    sData.h[1] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[0]) * (1.0 - xi[1])
    sData.h[2] = -0.25 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[0] + xi[1])
    sData.h[3] = 0.5 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[1])
    sData.h[4] = -0.25 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[0] - xi[1])
    sData.h[5] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[0]) * (1.0 + xi[1])
    sData.h[6] = -0.25 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[0] - xi[1])
    sData.h[7] = 0.5 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[1])

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -0.25 * (-1.0 + xi[1]) * (2.0 * xi[0] + xi[1])
    sData.dhdxi[1, 0] = xi[0] * (-1.0 + xi[1])
    sData.dhdxi[2, 0] = 0.25 * (-1.0 + xi[1]) * (-2.0 * xi[0] + xi[1])
    sData.dhdxi[3, 0] = -0.5 * (1.0 + xi[1]) * (-1.0 + xi[1])
    sData.dhdxi[4, 0] = 0.25 * (1.0 + xi[1]) * (2.0 * xi[0] + xi[1])
    sData.dhdxi[5, 0] = -xi[0] * (1.0 + xi[1])
    sData.dhdxi[6, 0] = -0.25 * (1.0 + xi[1]) * (-2.0 * xi[0] + xi[1])
    sData.dhdxi[7, 0] = 0.5 * (1.0 + xi[1]) * (-1.0 + xi[1])

    sData.dhdxi[0, 1] = -0.25 * (-1.0 + xi[0]) * (xi[0] + 2.0 * xi[1])
    sData.dhdxi[1, 1] = 0.5 * (1.0 + xi[0]) * (-1.0 + xi[0])
    sData.dhdxi[2, 1] = 0.25 * (1.0 + xi[0]) * (-xi[0] + 2.0 * xi[1])
    sData.dhdxi[3, 1] = -xi[1] * (1.0 + xi[0])
    sData.dhdxi[4, 1] = 0.25 * (1.0 + xi[0]) * (xi[0] + 2.0 * xi[1])
    sData.dhdxi[5, 1] = -0.5 * (1.0 + xi[0]) * (-1.0 + xi[0])
    sData.dhdxi[6, 1] = -0.25 * (-1.0 + xi[0]) * (-xi[0] + 2.0 * xi[1])
    sData.dhdxi[7, 1] = xi[1] * (-1.0 + xi[0])

    return sData

#-------------------------------------

def getShapeQuad9(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 2D quadrilateral element with 9 nodes (Quad9).
    
    Computes shape functions and derivatives at a single integration point
    using tensor product of Line3 functions.
    
    Args:
        xi: Location of integration point as ndarray of length 2
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 2D (length != 2)
    """

    if len(xi) != 2:
        raise NotImplementedError('2D only')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(9)
    sData.dhdxi = np.empty(shape=(9, 2))
    sData.xi = xi

    nodeMap = np.array([[0, 1, 2], [7, 8, 3], [6, 5, 4]])

    s0 = getShapeLine3(xi[0])
    s1 = getShapeLine3(xi[1])

    for i in range(3):
        for j in range(3):
            iNod = nodeMap[i, j]

            sData.h[iNod] = s0.h[i] * s1.h[j]
            sData.dhdxi[iNod, 0] = s0.h[i] * s1.dhdxi[0, j]
            sData.dhdxi[iNod, 1] = s0.dhdxi[0, i] * s1.h[j]

    return sData

#----------------------------------------------------------------------

def getShapeTetra4(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 3D tetrahedral element with 4 nodes (Tetra4).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of integration point as ndarray of length 3
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 3D (length != 3)
    """

    if len(xi) != 3:
        raise NotImplementedError('3D only')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(4)
    sData.dhdxi = np.empty(shape=(4, 3))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 1.0 - xi[0] - xi[1] - xi[2]
    sData.h[1] = xi[0]
    sData.h[2] = xi[1]
    sData.h[3] = xi[2]

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -1.0
    sData.dhdxi[1, 0] = 1.0
    sData.dhdxi[2, 0] = 0.0
    sData.dhdxi[3, 0] = 0.0

    sData.dhdxi[0, 1] = -1.0
    sData.dhdxi[1, 1] = 0.0
    sData.dhdxi[2, 1] = 1.0
    sData.dhdxi[3, 1] = 0.0

    sData.dhdxi[0, 2] = -1.0
    sData.dhdxi[1, 2] = 0.0
    sData.dhdxi[2, 2] = 0.0
    sData.dhdxi[3, 2] = 1.0

    return sData

#----------------------------------------------------------------------

def getShapePyramid5(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 3D pyramid element with 5 nodes (Pyramid5).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of integration point as ndarray of length 3
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 3D (length != 3)
    """

    if len(xi) != 3:
        raise NotImplementedError('3D only')

    sData = shapeData()

    # Set length of lists
    sData.h = np.empty(5)
    sData.dhdxi = np.empty(shape=(5, 3))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.h[1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.h[2] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.h[3] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.h[4] = 0.5 * (1.0 + xi[2])

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.dhdxi[1, 0] = 0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.dhdxi[2, 0] = 0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.dhdxi[3, 0] = -0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.dhdxi[4, 0] = 0.0

    sData.dhdxi[0, 1] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    sData.dhdxi[1, 1] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    sData.dhdxi[2, 1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    sData.dhdxi[3, 1] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    sData.dhdxi[4, 1] = 0.0

    sData.dhdxi[0, 2] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[1])
    sData.dhdxi[1, 2] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[1])
    sData.dhdxi[2, 2] = -0.125 * (1.0 + xi[0]) * (1.0 + xi[1])
    sData.dhdxi[3, 2] = -0.125 * (1.0 - xi[0]) * (1.0 + xi[1])
    sData.dhdxi[4, 2] = 0.5

    return sData

#----------------------------------------------------------------------

def getShapePrism6(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 3D prismatic element with 6 nodes (Prism6).
    
    Computes shape functions using tensor product of Line2 and Tria3 functions.
    
    Args:
        xi: Location of integration point as ndarray of length 3
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 3D (length != 3)
    """

    if len(xi) != 3:
        raise NotImplementedError('3D only')

    sData = shapeData()

    sData.h = np.empty(6)
    sData.dhdxi = np.empty(shape=(6, 3))
    sData.xi = xi

    sDataLine2 = getShapeLine2(xi[2])
    sDataTria3 = getShapeTria3(xi[:2])

    for i in range(3):
        for j in range(2):
            sData.h[i * 2 + j] = sDataLine2.h[j] * sDataTria3.h[i]

            sData.dhdxi[i * 2 + j, 0] = sDataLine2.h[j] * sDataTria3.dhdxi[i, 0]
            sData.dhdxi[i * 2 + j, 1] = sDataLine2.h[j] * sDataTria3.dhdxi[i, 1]
            sData.dhdxi[i * 2 + j, 2] = sDataLine2.dhdxi[j, 0] * sDataTria3.h[i]

    return sData

#----------------------------------------------------------------------

def getShapePrism18(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 3D prismatic element with 18 nodes (Prism18).
    
    Computes shape functions using tensor product of Line3 and Tria6 functions.
    
    Args:
        xi: Location of integration point as ndarray of length 3
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 3D (length != 3)
    """

    if len(xi) != 3:
        raise NotImplementedError('3D only')

    sData = shapeData()

    sData.h = np.empty(18)
    sData.dhdxi = np.empty(shape=(18, 3))
    sData.xi = xi

    sDataLine3 = getShapeLine3(xi[2])
    sDataTria6 = getShapeTria6(xi[:2])

    for i in range(6):
        for j in range(3):
            sData.h[i * 3 + j] = sDataLine3.h[j] * sDataTria6.h[i]

            sData.dhdxi[i * 3 + j, 0] = sDataLine3.h[j] * sDataTria6.dhdxi[i, 0]
            sData.dhdxi[i * 3 + j, 1] = sDataLine3.h[j] * sDataTria6.dhdxi[i, 1]
            sData.dhdxi[i * 3 + j, 2] = sDataLine3.dhdxi[0,j] * sDataTria6.h[i]

    return sData

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def getShapeHexa8(xi: np.ndarray) -> shapeData:
    """
    Shape functions for 3D hexahedral element with 8 nodes (Hexa8).
    
    Computes shape functions and derivatives at a single integration point
    in the parent coordinate system.
    
    Args:
        xi: Location of integration point as ndarray of length 3
        
    Returns:
        shapeData: Object containing h, dhdxi, and xi
        
    Raises:
        NotImplementedError: If input is not 3D (length != 3)
    """

    if len(xi) != 3:
        raise NotImplementedError('The isoparametric coordinate should be 3D.')

    sData = shapeData()

    sData.h = np.empty(8)
    sData.dhdxi = np.empty(shape=(8, 3))
    sData.xi = xi

    # Calculate shape functions
    sData.h[0] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.h[1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.h[2] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.h[3] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.h[4] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
    sData.h[5] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2])
    sData.h[6] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])
    sData.h[7] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2])

    # Calculate derivatives of shape functions
    sData.dhdxi[0, 0] = -0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.dhdxi[1, 0] = 0.125 * (1.0 - xi[1]) * (1.0 - xi[2])
    sData.dhdxi[2, 0] = 0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.dhdxi[3, 0] = -0.125 * (1.0 + xi[1]) * (1.0 - xi[2])
    sData.dhdxi[4, 0] = -0.125 * (1.0 - xi[1]) * (1.0 + xi[2])
    sData.dhdxi[5, 0] = 0.125 * (1.0 - xi[1]) * (1.0 + xi[2])
    sData.dhdxi[6, 0] = 0.125 * (1.0 + xi[1]) * (1.0 + xi[2])
    sData.dhdxi[7, 0] = -0.125 * (1.0 + xi[1]) * (1.0 + xi[2])

    sData.dhdxi[0, 1] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    sData.dhdxi[1, 1] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    sData.dhdxi[2, 1] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[2])
    sData.dhdxi[3, 1] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[2])
    sData.dhdxi[4, 1] = -0.125 * (1.0 - xi[0]) * (1.0 + xi[2])
    sData.dhdxi[5, 1] = -0.125 * (1.0 + xi[0]) * (1.0 + xi[2])
    sData.dhdxi[6, 1] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[2])
    sData.dhdxi[7, 1] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[2])

    sData.dhdxi[0, 2] = -0.125 * (1.0 - xi[0]) * (1.0 - xi[1])
    sData.dhdxi[1, 2] = -0.125 * (1.0 + xi[0]) * (1.0 - xi[1])
    sData.dhdxi[2, 2] = -0.125 * (1.0 + xi[0]) * (1.0 + xi[1])
    sData.dhdxi[3, 2] = -0.125 * (1.0 - xi[0]) * (1.0 + xi[1])
    sData.dhdxi[4, 2] = 0.125 * (1.0 - xi[0]) * (1.0 - xi[1])
    sData.dhdxi[5, 2] = 0.125 * (1.0 + xi[0]) * (1.0 - xi[1])
    sData.dhdxi[6, 2] = 0.125 * (1.0 + xi[0]) * (1.0 + xi[1])
    sData.dhdxi[7, 2] = 0.125 * (1.0 - xi[0]) * (1.0 + xi[1])

    return sData

#----------------------------------------------------------------------

def getElemType(elemCoords: np.ndarray) -> str:
    """
    Determine element type from nodal coordinates.
    
    Returns the element type name based on the shape and number of nodes
    in the element.
    
    Args:
        elemCoords: Array of nodal coordinates with shape (nNodes, nDimensions)
        
    Returns:
        str: Element type name
            - 1D: 'Line2', 'Line3'
            - 2D: 'Tria3', 'Tria6', 'Quad4', 'Quad8', 'Quad9'
            - 3D: 'Tetra4', 'Pyramid5', 'Prism6', 'Prism18', 'Hexa8'
            
    Raises:
        NotImplementedError: If element type not supported or rank is not 1, 2, or 3
    """

    nNel = elemCoords.shape[0]
    rank = elemCoords.shape[1]

    if rank == 1:
        if nNel == 2:
            return "Line2"
        elif nNel == 3:
            return "Line3"
        else:
            raise NotImplementedError('No 1D element with ' + str(nNel) + ' nodes available')
    elif rank == 2:
        if nNel == 3:
            return "Tria3"
        elif nNel == 4:
            return "Quad4"
        elif nNel == 6:
            return "Tria6"
        elif nNel == 8:
            return "Quad8"
        elif nNel == 9:
            return "Quad9"
        else:
            raise NotImplementedError('No 2D element with ' + str(nNel) + ' nodes available')
    elif rank == 3:
        if nNel == 4:
            return "Tetra4"
        elif nNel == 5:
            return "Pyramid5"
        elif nNel == 6:
            return "Prism6"
        elif nNel == 8:
            return "Hexa8"
        elif nNel == 18:
            return "Prism18"
        else:
            raise NotImplementedError('No 3D element with ' + str(nNel) + ' nodes available')
    else:
        raise NotImplementedError('Rank must be 1, 2, or 3')

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def tria_scheme(order: int) -> Tuple[List[List[float]], List[float]]:
    """
    Integration scheme for 2D triangular elements.
    
    Returns integration point coordinates and weights for triangular elements.
    
    Args:
        order: Integration order (1, 3, or 7 integration points)
        
    Returns:
        Tuple of:
            - List of integration point coordinates in parent element
            - List of corresponding integration weights
            
    Raises:
        NotImplementedError: If order is not 1, 3, or 7
    """

    if order == 1:
        xi = [[1.0 / 3.0, 1.0 / 3.0]]
        weight = [0.5]
    elif order == 3:
        r1 = 1.0 / 6.0
        r2 = 2.0 / 3.0

        xi = [[r1, r1], [r2, r1], [r1, r2]]

        w1 = 1.0 / 6.0

        weight = [w1, w1, w1]
    elif order == 7:
        r1 = 0.5 * 0.1012865073235
        r2 = 0.5 * 0.7974269853531
        r4 = 0.5 * 0.4701420641051
        r6 = 0.0597158717898
        r7 = 1.0 / 3.0

        xi = [[r1, r1], [r2, r1], [r1, r2], [r4, r6], [r4, r4], [r6, r4], [r7, r7]]

        w1 = 0.1259391805448
        w4 = 0.1323941527885
        w7 = 0.225

        weight = [w1, w1, w1, w4, w4, w4, w7]
    else:
        raise NotImplementedError('Order must be 1, 3, or 7')

    return xi, weight

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def tetra_scheme(order: int) -> Tuple[List[List[float]], List[float]]:
    """
    Integration scheme for 3D tetrahedral elements.
    
    Returns integration point coordinates and weights for tetrahedral elements.
    
    Args:
        order: Integration order (1 for single integration point)
        
    Returns:
        Tuple of:
            - List of integration point coordinates in parent element
            - List of corresponding integration weights
            
    Raises:
        NotImplementedError: If order is not 1
    """

    if order == 1:
        third = 1.0 / 3.0

        xi = [[third, third, third]]
        weight = [0.5 * third]
    else:
        raise NotImplementedError('Only order 1 integration implemented')

    return xi, weight

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def pyramid_scheme(order: int) -> Tuple[List[List[float]], List[float]]:
    """
    Integration scheme for 3D pyramid elements.
    
    Returns integration point coordinates and weights for pyramid elements.
    
    Args:
        order: Integration order (1 for single integration point)
        
    Returns:
        Tuple of:
            - List of integration point coordinates in parent element
            - List of corresponding integration weights
            
    Raises:
        NotImplementedError: If order is not 1
    """

    if order == 1:
        xi = [[0., 0., -0.5]]
        weight = [128.0 / 27.0]
    else:
        raise NotImplementedError('Only order 1 integration implemented')

    return xi, weight

#-----------------------------------------------------------------------

def getIntegrationPoints(elemType: str, order: int, scheme: str) -> Tuple[List[List[float]], List[float]]:
    """
    Get integration points and weights for any element type.
    
    Returns integration scheme (coordinates and weights) for any supported element
    type. Allows modification of standard integration order.
    
    Args:
        elemType: Element type name
            - 1D: 'Line2', 'Line3'
            - 2D: 'Tria3', 'Tria6', 'Quad4', 'Quad8', 'Quad9'
            - 3D: 'Tetra4', 'Pyramid5', 'Prism6', 'Prism18', 'Hexa8'
        order: Integration order adjustment
            - 0: Standard integration for element type
            - +1: Higher order integration
            - -1: Lower order integration
        scheme: Integration scheme name (e.g., 'Gauss')
        
    Returns:
        Tuple of:
            - List of integration point coordinates
            - List of integration weights
            
    Raises:
        NotImplementedError: If element type is unknown
    """

    xi = []
    weight = []

    if elemType[:-1] == "Line":
        if elemType == "Line2":
            stdOrder = 2
        elif elemType == "Line3":
            stdOrder = 3
        xi, weight = gauss_scheme(stdOrder + order)
        xi = [float(a.real) for a in xi]

    elif elemType[:-1] == "Tria":
        orderArray = [1, 3, 7]
        if elemType == "Tria3":
            stdOrder = 0
        elif elemType == "Tria6":
            stdOrder = 1
        xi, weight = tria_scheme(orderArray[stdOrder + order])

    elif elemType[:-1] == "Tetra":
        stdOrder = 1
        xi, weight = tetra_scheme(stdOrder + order)

    elif elemType == "Pyramid5":
        stdOrder = 1
        xi, weight = pyramid_scheme(stdOrder + order)

    elif elemType[:-1] == "Quad":
        if elemType == "Quad4":
            stdOrder = 2
        elif elemType == "Quad8" or elemType == "Quad9":
            stdOrder = 3
        stdOrder += order

        ip, w = gauss_scheme(stdOrder)

        for i in range(stdOrder):
            for j in range(stdOrder):
                xi.append([float(ip[i].real), float(ip[j].real)])
                weight.append(w[i] * w[j])

    elif elemType[:-1] == "Hexa":
        if elemType == "Hexa8":
            stdOrder = 2

        stdOrder += order

        ip, w = gauss_scheme(stdOrder)

        for i in range(stdOrder):
            for j in range(stdOrder):
                for k in range(stdOrder):
                    xi.append([float(ip[i].real), float(ip[j].real), float(ip[k].real)])
                    weight.append(w[i] * w[j] * w[k])

    elif elemType[:-1] == "Prism":
        orderArray = [1, 3, 7]
        if elemType == "Prism6":
            stdOrder = 2
        elif elemType == "Prism18":
            stdOrder = 3

        stdOrder += order

        ip0, w0 = tria_scheme(orderArray[stdOrder])
        ip1, w1 = gauss_scheme(stdOrder)

        for i in range(orderArray[stdOrder]):
            for j in range(stdOrder):
                xi.append([float(ip0[i][0]), float(ip0[i][1]), float(ip1[j].real)])
                weight.append(w0[i] * w1[j])
    else:
        raise NotImplementedError('Element type not known.')

    return xi, weight

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def calcWeightandDerivatives(elemCoords: np.ndarray, sData: shapeData, weight: float) -> None:
    """
    Calculate physical derivatives and integration weight.
    
    Computes the derivatives of shape functions with respect to physical coordinates
    and the weighted integration weight using the Jacobian transformation.
    Results are stored in sData.
    
    Args:
        elemCoords: Array of nodal coordinates (nNodes, nDimensions)
        sData: Shape data object containing dhdxi and xi
        weight: Integration weight in parent element
        
    Modifies:
        sData.dhdx: Derivatives w.r.t. physical coordinates
        sData.weight: Integration weight in physical element
    """

    jac = np.dot(elemCoords.transpose(), sData.dhdxi)

    if jac.shape[0] == jac.shape[1]:
        sData.dhdx = np.dot(sData.dhdxi, inv(jac))
        sData.weight = abs(det(jac)) * weight

    elif jac.shape[0] == 2 and jac.shape[1] == 1:
        sData.weight = sqrt(sum(sum(jac * jac))) * weight

    elif jac.shape[0] == 3 and jac.shape[1] == 2:
        jac3 = np.zeros(shape=(3, 3))

        jac3[:, :2] = jac

        dA = np.zeros(3)

        dA[0] = norm(np.cross(jac3[:, 1], jac3[:, 2]))
        dA[1] = norm(np.cross(jac3[:, 2], jac3[:, 0]))
        dA[2] = norm(np.cross(jac3[:, 0], jac3[:, 1]))

        sData.weight = norm(dA) * weight

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def getElemShapeData(elemCoords: np.ndarray, order: int = 0,
                     method: str = 'Gauss', elemType: str = 'Default') -> elemShapeData:
    """
    Get shape function and integration data for an element.
    
    Determines element type, generates integration points, and computes
    shape functions and physical derivatives at each integration point.
    
    Args:
        elemCoords: Array of nodal coordinates (nNodes, nDimensions)
        order: Integration order adjustment (default 0)
        method: Integration scheme name, e.g., 'Gauss' (default)
        elemType: Element type name. If 'Default', determined from elemCoords
        
    Returns:
        elemShapeData: Object containing shape data for all integration points
        
    Raises:
        NotImplementedError: If element type is unknown
    """

    elemData = elemShapeData()

    if elemType == 'Default':
        elemType = getElemType(elemCoords)

    (intCrds, intWghts) = getIntegrationPoints(elemType, order, method)

    for xi, weight in zip(intCrds, intWghts):
        try:
            sData = eval('getShape' + elemType + '(xi)')
        except Exception:
            raise NotImplementedError('Unknown type: ' + elemType)

        calcWeightandDerivatives(elemCoords, sData, weight)

        sData.x = np.dot(sData.h, elemCoords)

        elemData.sData.append(sData)

    return elemData

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------

def getShapeData(order: int = 0, method: str = 'Gauss', elemType: str = 'Default') -> elemShapeData:
    """
    Get shape function and integration data for an element type.
    
    Generates integration points and computes shape functions and parent element
    derivatives at each integration point. Does not compute physical derivatives.
    
    Args:
        order: Integration order adjustment (default 0)
        method: Integration scheme name, e.g., 'Gauss' (default)
        elemType: Element type name (required)
        
    Returns:
        elemShapeData: Object containing shape data for all integration points
        
    Raises:
        NotImplementedError: If element type is unknown
    """

    shpData = elemShapeData()

    (intCrds, intWghts) = getIntegrationPoints(elemType, order, method)

    for xi, weight in zip(intCrds, intWghts):
        try:
            sData = eval('getShape' + elemType + '(xi)')
        except Exception:
            raise NotImplementedError('Unknown type: ' + elemType)

        sData.dhdx = sData.dhdxi
        sData.weight = weight

        shpData.sData.append(sData)

    return shpData
