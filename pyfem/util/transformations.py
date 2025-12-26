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

from numpy import array, dot, ndarray, empty, ix_
from scipy.linalg import norm
from typing import Union

from .dataStructures import GlobalData


def getRotationMatrix(el_coords: ndarray) -> ndarray:
    """
    Compute the 2D rotation matrix for an element.
    
    Calculates the rotation matrix that transforms from global coordinates to
    element coordinates based on the element's orientation in 2D space.
    
    Args:
        el_coords: Element node coordinates as (2, 2) array where rows are nodes
        
    Returns:
        2x2 rotation matrix for coordinate transformation
        
    Raises:
        NotImplementedError: If element coordinates are not 2D (shape[1] != 2)
    """

    # Check the dimension of physical space
    if el_coords.shape[1] != 2:
        raise NotImplementedError('Rotation matrix only implemented for 2D situation')

    # Compute the (undeformed) element length
    l0 = norm(el_coords[1] - el_coords[0])

    # Set up the rotation matrix to rotate a global
    # coordinate to an element coordinate (see Ch 1.3)
    sinalpha = (el_coords[1, 1] - el_coords[0, 1]) / l0
    cosalpha = (el_coords[1, 0] - el_coords[0, 0]) / l0

    return array([[cosalpha, sinalpha], [-sinalpha, cosalpha]])


def vectorToElementCoordinates(a: ndarray, el_coords: ndarray) -> ndarray:
    """
    Transform a vector from global to element coordinates.
    
    Rotates a vector (or block of vectors) from the global coordinate system
    to the element coordinate system. Handles vectors with multiple components.
    
    Args:
        a: Vector in global coordinates with shape (n,) where n is divisible by 2
        el_coords: Element node coordinates as (2, 2) array
        
    Returns:
        Vector in element coordinates with same shape as input
        
    Raises:
        RuntimeError: If vector length is not divisible by rotation matrix size
    """

    R = getRotationMatrix(el_coords)

    a_bar = empty(a.shape)

    if len(a_bar) % len(R) != 0:
        raise RuntimeError('Vector does not have the right shape to be rotated')

    for i in range(len(a_bar) // len(R)):
        a_bar[len(R) * i:len(R) * (i + 1)] = dot(R, a[len(R) * i:len(R) * (i + 1)])

    return a_bar


def matrixToElementCoordinates(a: ndarray, el_coords: ndarray) -> ndarray:
    """
    Transform a matrix from global to element coordinates.
    
    Rotates a square matrix (or block matrix) from the global coordinate system
    to the element coordinate system using: a_bar = R @ a @ R^T
    
    Args:
        a: Square matrix in global coordinates with shape (m, m) 
           where m is divisible by 2
        el_coords: Element node coordinates as (2, 2) array
        
    Returns:
        Matrix in element coordinates with same shape as input
        
    Raises:
        RuntimeError: If matrix dimensions are not divisible by rotation matrix size
    """

    R = getRotationMatrix(el_coords)

    a_bar = empty(a.shape)

    if a_bar.shape[0] % len(R) != 0 or a_bar.shape[1] % len(R) != 0:
        raise RuntimeError('Matrix does not have the right shape to be rotated')

    for i in range(a_bar.shape[0] // len(R)):
        iran = list(range(len(R) * i, len(R) * (i + 1)))

        for j in range(a_bar.shape[1] // len(R)):
            jran = list(range(len(R) * j, len(R) * (j + 1)))

            a_bar[ix_(iran, jran)] = dot(dot(R, a[ix_(iran, jran)]), R.transpose())

    return a_bar


def vectorToGlobalCoordinates(a_bar: ndarray, el_coords: ndarray) -> ndarray:
    """
    Transform a vector from element to global coordinates.
    
    Rotates a vector (or block of vectors) from the element coordinate system
    back to the global coordinate system. Inverse of vectorToElementCoordinates.
    
    Args:
        a_bar: Vector in element coordinates with shape (n,) where n is divisible by 2
        el_coords: Element node coordinates as (2, 2) array
        
    Returns:
        Vector in global coordinates with same shape as input
        
    Raises:
        RuntimeError: If vector length is not divisible by rotation matrix size
    """

    R = getRotationMatrix(el_coords)

    a = empty(a_bar.shape)

    if len(a) % len(R) != 0:
        raise RuntimeError('Vector does not have the right shape to be rotated')

    for i in range(len(a) // len(R)):
        a[len(R) * i:len(R) * (i + 1)] = dot(R.transpose(), a_bar[len(R) * i:len(R) * (i + 1)])

    return a


def matrixToGlobalCoordinates(a_bar: ndarray, el_coords: ndarray) -> ndarray:
    """
    Transform a matrix from element to global coordinates.
    
    Rotates a square matrix (or block matrix) from the element coordinate system
    back to the global coordinate system using: a = R^T @ a_bar @ R
    Inverse of matrixToElementCoordinates.
    
    Args:
        a_bar: Square matrix in element coordinates with shape (m, m)
               where m is divisible by 2
        el_coords: Element node coordinates as (2, 2) array
        
    Returns:
        Matrix in global coordinates with same shape as input
        
    Raises:
        RuntimeError: If matrix dimensions are not divisible by rotation matrix size
    """

    R = getRotationMatrix(el_coords)

    a = empty(a_bar.shape)

    if a.shape[0] % len(R) != 0 or a.shape[1] % len(R) != 0:
        raise RuntimeError('Matrix does not have the right shape to be rotated')

    for i in range(a.shape[0] // len(R)):
        iran = list(range(len(R) * i, len(R) * (i + 1)))

        for j in range(a.shape[1] // len(R)):
            jran = list(range(len(R) * j, len(R) * (j + 1)))

            a[ix_(iran, jran)] = dot(dot(R.transpose(), a_bar[ix_(iran, jran)]), R)

    return a


def toElementCoordinates(a: Union[ndarray], el_coords: ndarray) -> ndarray:
    """
    Transform array (vector or matrix) from global to element coordinates.
    
    Dispatcher function that automatically detects input type (1D vector or 2D matrix)
    and applies the appropriate transformation. Rotates from global to element 
    coordinate system.
    
    Args:
        a: Input vector (shape (n,)) or matrix (shape (m, m)) in global coordinates
        el_coords: Element node coordinates as (2, 2) array
        
    Returns:
        Transformed array in element coordinates with same shape as input
        
    Raises:
        NotImplementedError: If input is not a 1D or 2D ndarray
        RuntimeError: If array dimensions are incompatible with rotation matrix
    """

    # Vector
    if isinstance(a, ndarray) and len(a.shape) == 1:
        return vectorToElementCoordinates(a, el_coords)

    # Matrix
    elif isinstance(a, ndarray) and len(a.shape) == 2:
        return matrixToElementCoordinates(a, el_coords)

    # Error
    else:
        raise NotImplementedError('Rotation to element coordinate system only works for vectors and matrices.')


def toGlobalCoordinates(a: Union[ndarray], el_coords: ndarray) -> ndarray:
    """
    Transform array (vector or matrix) from element to global coordinates.
    
    Dispatcher function that automatically detects input type (1D vector or 2D matrix)
    and applies the appropriate transformation. Rotates from element to global
    coordinate system. Inverse of toElementCoordinates.
    
    Args:
        a: Input vector (shape (n,)) or matrix (shape (m, m)) in element coordinates
        el_coords: Element node coordinates as (2, 2) array
        
    Returns:
        Transformed array in global coordinates with same shape as input
        
    Raises:
        NotImplementedError: If input is not a 1D or 2D ndarray
        RuntimeError: If array dimensions are incompatible with rotation matrix
    """

    # Vector
    if isinstance(a, ndarray) and len(a.shape) == 1:
        return vectorToGlobalCoordinates(a, el_coords)

    # Matrix
    elif isinstance(a, ndarray) and len(a.shape) == 2:
        return matrixToGlobalCoordinates(a, el_coords)

    # Error
    else:
        raise NotImplementedError('Rotation to global coordinate system only works for vectors and matrices.')
