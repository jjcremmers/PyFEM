# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

from typing import Tuple
from pyfem.materials.BaseMaterial import BaseMaterial
import numpy as np


class Isotropic(BaseMaterial):
    """
    Isotropic linear elastic material model.
    
    This class implements a standard isotropic linear elastic material model
    using Hooke's law. The material behavior is fully characterized by two
    elastic constants: Young's modulus (E) and Poisson's ratio (nu).
    
    The model can operate in two modes:
    - Total formulation: stress = H · strain
    - Incremental formulation: Δσ = H · Δε
    
    Attributes
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.
    H : ndarray
        6x6 elastic stiffness matrix (Hookean matrix).
    incremental : bool
        Flag for incremental or total formulation.
    outLabels : list of str
        Labels for stress components.
    
    Notes
    -----
    The stress and strain vectors follow Voigt notation:
    [σ11, σ22, σ33, σ23, σ13, σ12]
    """

    def __init__(self, props) -> None:
        """
        Initialize the isotropic material model.
        
        Constructs the 6x6 elastic stiffness matrix based on Young's modulus
        and Poisson's ratio.
        
        Parameters
        ----------
        props : object
            Properties object containing material parameters. Must have:
            - E : float
                Young's modulus.
            - nu : float
                Poisson's ratio.
            - incremental : bool, optional
                Flag for incremental formulation (default: False).
        """
        self.incremental = False
        
        # Call the BaseMaterial constructor
        BaseMaterial.__init__(self, props)

        # Create the hookean matrix
        self.H = np.zeros((6, 6))

        fac = 1.0 / (2.0 * self.nu * self.nu + self.nu - 1.0)
      
        self.H[0, 0] = fac * self.E * (self.nu - 1.0)
        self.H[0, 1] = -1.0 * fac * self.E * self.nu
        self.H[0, 2] = self.H[0, 1]
        self.H[1, 0] = self.H[0, 1]
        self.H[1, 1] = self.H[0, 0]
        self.H[1, 2] = self.H[0, 1]
        self.H[2, 0] = self.H[0, 1]
        self.H[2, 1] = self.H[0, 1]
        self.H[2, 2] = self.H[0, 0]
        self.H[3, 3] = self.E / (2.0 + 2.0 * self.nu)
        self.H[4, 4] = self.H[3, 3]
        self.H[5, 5] = self.H[3, 3]

        # Set the labels for the output data in this material model
        self.outLabels = ["S11", "S22", "S33", "S23", "S13", "S12"]
        
        if self.incremental:
            self.setHistoryParameter('sigma', np.zeros(6))
            self.commitHistory()

    # ---------------------------------------------------------------------------
    #
    # ---------------------------------------------------------------------------

    def getStress(self, deformation) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute stress and material tangent matrix.
        
        Parameters
        ----------
        deformation : object
            Deformation object containing strain information. Must have:
            - strain : ndarray
                Total strain vector (for total formulation).
            - dstrain : ndarray
                Strain increment vector (for incremental formulation).
        
        Returns
        -------
        sigma : ndarray
            Stress vector in Voigt notation [σ11, σ22, σ33, σ23, σ13, σ12].
        H : ndarray
            Material tangent stiffness matrix (6x6).
        """
        if self.incremental:
            sigma = self.getHistoryParameter('sigma')
            sigma += self.H @ deformation.dstrain
            self.setHistoryParameter('sigma', sigma)
        else:
            sigma = self.H @ deformation.strain

        self.outData = sigma

        return sigma, self.H

    def getTangent(self) -> np.ndarray:
        """
        Get the material tangent stiffness matrix.
        
        Returns
        -------
        ndarray
            Material tangent stiffness matrix H (6x6).
        """
        return self.H

