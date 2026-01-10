# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""
Tsai-Wu failure criterion for composite materials.

This module implements the Tsai-Wu failure criterion, which is a tensor
polynomial criterion that accounts for interaction between stress components.
It is widely used for orthotropic materials such as fiber-reinforced composites.
"""

from pyfem.materials.BaseFailure import BaseFailure
from numpy import ndarray, sqrt
from typing import Union


class TsaiWu(BaseFailure):
    """
    Tsai-Wu failure criterion for orthotropic materials.
    
    The Tsai-Wu criterion is a general tensor polynomial failure theory that
    accounts for interaction between stress components. It is expressed as:
    
    FI = F₁σ₁ + F₂σ₂ + F₆τ₁₂ + F₁₁σ₁² + F₂₂σ₂² + F₆₆τ₁₂² + 2F₁₂σ₁σ₂
    
    where the coefficients are defined in terms of strength values:
    - F₁ = 1/Xᵗ - 1/Xᶜ
    - F₂ = 1/Yᵗ - 1/Yᶜ
    - F₁₁ = 1/(Xᵗ·Xᶜ)
    - F₂₂ = 1/(Yᵗ·Yᶜ)
    - F₆₆ = 1/S²
    - F₁₂ = interaction coefficient (typically -0.5·sqrt(F₁₁·F₂₂))
    
    Failure occurs when FI ≥ 1.0
    
    Required Properties
    -------------------
    Xt : float
        Tensile strength in fiber direction (direction 1).
    Xc : float
        Compressive strength in fiber direction (positive value).
    Yt : float
        Tensile strength transverse to fibers (direction 2).
    Yc : float
        Compressive strength transverse to fibers (positive value).
    S : float
        Shear strength.
    F12 : float, optional
        Interaction coefficient F₁₂. If not provided, computed as:
        F₁₂ = -0.5 × sqrt(F₁₁ × F₂₂)
        This is a typical assumption when experimental data is unavailable.
    
    Examples
    --------
    Typical values for carbon fiber/epoxy composite:
    
    <Failure>
        type = "TsaiWu";
        Xt = 1500.0;    # MPa - fiber tensile strength
        Xc = 1200.0;    # MPa - fiber compressive strength
        Yt = 50.0;      # MPa - transverse tensile strength
        Yc = 200.0;     # MPa - transverse compressive strength
        S = 100.0;      # MPa - shear strength
        F12 = -0.5e-6;  # Optional interaction coefficient (1/MPa²)
    </Failure>
    
    Notes
    -----
    Advantages:
    - Accounts for interaction between stress components
    - Single equation for all loading conditions
    - Distinguishes between tension and compression
    - Widely validated for composite materials
    
    Limitations:
    - Requires interaction coefficient F₁₂ (often unknown)
    - Does not predict failure mode
    - May be non-conservative for some loading cases
    
    The criterion reduces to von Mises for isotropic materials with appropriate
    strength values.
    """

    def __init__(self, props) -> None:
        """
        Initialize the Tsai-Wu failure criterion.
        
        Parameters
        ----------
        props : Properties
            Material properties object containing strength values.
        """
        BaseFailure.__init__(self, props)
        
        # Validate required properties
        required_props = ['Xt', 'Xc', 'Yt', 'Yc', 'S']
        for prop in required_props:
            if not hasattr(self, prop):
                raise ValueError(f"Required property '{prop}' not found for TsaiWu failure criterion")

        # Compute Tsai-Wu coefficients
        self.F1 = 1.0 / self.Xt - 1.0 / self.Xc
        self.F2 = 1.0 / self.Yt - 1.0 / self.Yc
        self.F11 = 1.0 / (self.Xt * self.Xc)
        self.F22 = 1.0 / (self.Yt * self.Yc)
        self.F66 = 1.0 / (self.S * self.S)
        
        # Interaction coefficient
        # If not provided, use typical approximation: F12 = -0.5 * sqrt(F11 * F22)
        if not hasattr(self, 'F12'):
            self.F12 = -0.5 * sqrt(self.F11 * self.F22)
        
        # F6 coefficient (usually zero for orthotropic materials)
        self.F6 = 0.0

    def check(self, stress: ndarray, deformation) -> float:
        """
        Evaluate the Tsai-Wu failure criterion.
        
        Parameters
        ----------
        stress : ndarray
            Stress tensor in material principal directions.
            For 2D: [σ₁, σ₂, τ₁₂]
            For 3D: [σ₁, σ₂, σ₃, τ₂₃, τ₁₃, τ₁₂]
        deformation : object
            Deformation object (not used but kept for interface consistency).
        
        Returns
        -------
        float
            Failure index (FI). FI ≥ 1.0 indicates failure.
            The square root of FI gives the strength ratio.
        
        Notes
        -----
        The Tsai-Wu criterion is evaluated in the material coordinate system,
        so stress components must be rotated to principal material directions
        before calling this function.
        
        For plane stress problems (2D), only σ₁, σ₂, and τ₁₂ are considered.
        """
        
        # Extract stress components
        sigma1 = stress[0]
        sigma2 = stress[1]
        
        # Shear stress component
        if len(stress) == 3:
            tau12 = stress[2]  # 2D plane stress
        else:
            tau12 = stress[5]  # 3D (τ₁₂ is 6th component)

        # Compute Tsai-Wu failure index
        # FI = F₁σ₁ + F₂σ₂ + F₆τ₁₂ + F₁₁σ₁² + F₂₂σ₂² + F₆₆τ₁₂² + 2F₁₂σ₁σ₂
        FI = (self.F1 * sigma1 + 
              self.F2 * sigma2 + 
              self.F6 * tau12 +
              self.F11 * sigma1 * sigma1 +
              self.F22 * sigma2 * sigma2 +
              self.F66 * tau12 * tau12 +
              2.0 * self.F12 * sigma1 * sigma2)

        # Return square root for strength ratio interpretation
        # If FI > 0, return sqrt(FI), otherwise return -sqrt(|FI|)
        if FI >= 0.0:
            return sqrt(FI)
        else:
            # Negative FI can occur far from failure surface
            # Return negative value to indicate safe condition
            return -sqrt(abs(FI))
