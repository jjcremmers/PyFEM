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


class Kinematics:
    """
    Container for kinematic state variables at a material point.

    This class stores the kinematic state of a material point in a finite element
    analysis, including the deformation gradient tensor and various strain measures.
    The class is typically used within constitutive models to track the deformation
    and strain history of material points during nonlinear finite element analysis.

    Attributes:
        F (np.ndarray): Deformation gradient tensor of shape (nDim, nDim).
            Maps reference coordinates to spatial coordinates.
        E (np.ndarray): Green-Lagrange strain tensor of shape (nDim, nDim).
            Symmetric strain measure based on deformation gradient.
        strain (np.ndarray): Strain vector of shape (nStr,).
            Contains strain components in voigt notation (2, 3, or 6 components).
        dgdstrain (np.ndarray): Derivative of some scalar with respect to strain,
            shape (nStr,). Used in constitutive and material models.
        g (float): A scalar variable (e.g., equivalent strain, damage parameter).

    Args:
        nDim (int): Number of spatial dimensions of the problem. Must be 2 or 3.
        nStr (int): Number of strain components. Typically 2 (2D), 3 (axisymmetric),
            or 6 (3D full tensor in Voigt notation).

    Examples:
        Create a 3D kinematics object:
        >>> kin = Kinematics(nDim=3, nStr=6)
        >>> kin.F.shape
        (3, 3)
        >>> kin.strain.shape
        (6,)

        Create a 2D plane strain kinematics object:
        >>> kin = Kinematics(nDim=2, nStr=3)
        >>> kin.F.shape
        (2, 2)
    """

    def __init__(self, nDim: int, nStr: int) -> None:
        """
        Initialize the Kinematics object with zeros for all state variables.

        Args:
            nDim (int): Number of spatial dimensions (2 or 3).
            nStr (int): Number of strain components (2, 3, or 6).

        Raises:
            ValueError: If nDim is not 2 or 3, or if nStr is not valid.
        """
        if nDim not in (2, 3):
            raise ValueError(f"nDim must be 2 or 3, got {nDim}")
        if nStr not in (2, 3, 6):
            raise ValueError(f"nStr must be 2, 3, or 6, got {nStr}")

        self.F: np.ndarray = np.zeros(shape=(nDim, nDim))
        """Deformation gradient tensor."""

        self.E: np.ndarray = np.zeros(shape=(nDim, nDim))
        """Green-Lagrange strain tensor."""

        self.strain: np.ndarray = np.zeros(nStr)
        """Strain vector in Voigt notation."""

        self.dgdstrain: np.ndarray = np.zeros(nStr)
        """Derivative of scalar quantity with respect to strain."""

        self.g: float = 0.0
        """Scalar quantity (e.g., equivalent strain or damage parameter)."""
