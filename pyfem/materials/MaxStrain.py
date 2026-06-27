# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""
Maximum strain failure criterion for composite materials.

This module implements the maximum strain failure criterion, which is commonly
used for fiber-reinforced composite materials. Failure is predicted when any
principal strain component exceeds its allowable limit.
"""

from pyfem.materials.BaseFailure import BaseFailure
from numpy import ndarray
from typing import Union


class MaxStrain(BaseFailure):
    """
    Maximum strain failure criterion.
    
    This criterion predicts failure when any strain component exceeds its
    corresponding allowable strain. It is one of the simplest failure criteria
    for orthotropic materials like fiber-reinforced composites.
    
    The failure index is computed as:
    FI = max(ε₁/ε₁ᵗ, |ε₁|/ε₁ᶜ, ε₂/ε₂ᵗ, |ε₂|/ε₂ᶜ, |γ₁₂|/γ₁₂ᵘ)
    
    where:
    - ε₁, ε₂: normal strains in principal material directions
    - γ₁₂: shear strain
    - Superscripts: t = tension, c = compression, u = ultimate
    
    Failure occurs when FI ≥ 1.0
    
    Required Properties
    -------------------
    eps1t : float
        Maximum tensile strain in fiber direction (direction 1).
    eps1c : float
        Maximum compressive strain in fiber direction (magnitude, positive).
    eps2t : float
        Maximum tensile strain transverse to fibers (direction 2).
    eps2c : float
        Maximum compressive strain transverse to fibers (magnitude, positive).
    gamma12u : float
        Maximum shear strain (magnitude, positive).
    
    Examples
    --------
    Typical values for carbon fiber/epoxy composite:
    
    <Failure>
        type = "MaxStrain";
        eps1t = 0.015;      # 1.5% tensile strain in fiber direction
        eps1c = 0.012;      # 1.2% compressive strain
        eps2t = 0.005;      # 0.5% transverse tensile strain
        eps2c = 0.010;      # 1.0% transverse compressive strain
        gamma12u = 0.020;   # 2.0% shear strain
    </Failure>
    
    Notes
    -----
    - This criterion does not account for interaction between strain components
    - Conservative for combined loading conditions
    - Simple to implement and computationally efficient
    - Suitable for preliminary design and conservative estimates
    """

    def __init__(self, props) -> None:
        """
        Initialize the MaxStrain failure criterion.
        
        Parameters
        ----------
        props : Properties
            Material properties object containing strain limits.
        """
        BaseFailure.__init__(self, props)
        
        # Validate that all required properties are present
        required_props = ['eps1t', 'eps1c', 'eps2t', 'eps2c', 'gamma12u']
        for prop in required_props:
            if not hasattr(self, prop):
                raise ValueError(f"Required property '{prop}' not found for MaxStrain failure criterion")

    def check(self, stress: ndarray, deformation) -> float:
        """
        Evaluate the maximum strain failure criterion.
        
        Parameters
        ----------
        stress : ndarray
            Stress tensor (not used in this criterion, but kept for interface consistency).
        deformation : object
            Deformation object containing strain information.
            Must have attribute 'strain' or 'eps'.
        
        Returns
        -------
        float
            Failure index (FI). FI ≥ 1.0 indicates failure.
            FI < 1.0 indicates safe condition.
        
        Notes
        -----
        The strain tensor should be in engineering notation:
        - For 2D: [ε₁, ε₂, γ₁₂]
        - For 3D: [ε₁, ε₂, ε₃, γ₂₃, γ₁₃, γ₁₂]
        
        Shear strains are engineering shear strains (γ = 2ε₁₂).
        """
        # Get strain from deformation object
        if hasattr(deformation, 'strain'):
            strain = deformation.strain
        elif hasattr(deformation, 'eps'):
            strain = deformation.eps
        else:
            raise AttributeError("Deformation object must have 'strain' or 'eps' attribute")

        FI = 0.0

        # Check normal strain in direction 1 (fiber direction)
        if strain[0] > 0.0:
            FI = strain[0] / self.eps1t
        else:
            FI = abs(strain[0]) / self.eps1c

        # Check normal strain in direction 2 (transverse direction)
        if strain[1] > 0.0:
            FI = max(strain[1] / self.eps2t, FI)
        else:
            FI = max(abs(strain[1]) / self.eps2c, FI)

        # Check shear strain
        # For 2D plane stress: strain[2] is γ₁₂
        # For 3D: strain[5] is γ₁₂
        if len(strain) == 3:
            FI = max(abs(strain[2]) / self.gamma12u, FI)
        else:
            FI = max(abs(strain[5]) / self.gamma12u, FI)

        return FI
