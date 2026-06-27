# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Tuple
from pyfem.materials.BaseMaterial import BaseMaterial
import numpy as np


class Dummy(BaseMaterial):
    """
    Dummy material model for interface elements.
    
    This is a simple linear elastic material model used primarily for testing
    and interface elements. It provides a diagonal stiffness matrix scaled by
    a material parameter D.
    
    Attributes
    ----------
    H : ndarray
        Material stiffness matrix (diagonal).
    outLabels : list of str
        Labels for output variables (normal and tangential tractions).
    D : float
        Material stiffness parameter (inherited from props).
    
    Notes
    -----
    The material supports 2D (rank=2) and 3D (rank=3) interface elements with
    normal and tangential traction components.
    """

    def __init__(self, props) -> None:
        """
        Initialize the Dummy material model.
        
        Parameters
        ----------
        props : object
            Properties object containing material parameters. Must have:
            - rank : int
                Spatial dimension (2 or 3).
            - D : float
                Stiffness parameter.
        """
        BaseMaterial.__init__(self, props)

        if props.rank == 2:
            self.H = self.D * np.eye(2)
            self.outLabels = ["Tn", "Ts"]
        elif props.rank == 3:
            self.H = self.D * np.eye(3)
            self.outLabels = ["Tn", "Ts1", "Ts2"]

    def getStress(self, deformation) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute stress (traction) and material tangent matrix.
        
        Parameters
        ----------
        deformation : object
            Deformation object containing the strain field.
            Must have a strain attribute (ndarray).
        
        Returns
        -------
        sigma : ndarray
            Stress (traction) vector.
        H : ndarray
            Material tangent stiffness matrix.
        """
        sigma = self.H @ deformation.strain

        self.outData = sigma

        return sigma, self.H

    def getTangent(self) -> np.ndarray:
        """
        Get the material tangent stiffness matrix.
        
        Returns
        -------
        ndarray
            Material tangent stiffness matrix H.
        """
        return self.H

