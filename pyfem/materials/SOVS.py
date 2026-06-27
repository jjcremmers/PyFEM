# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""
Skorohod-Olevsky viscous sintering model for ceramic materials.

This module implements the Skorohod-Olevsky constitutive model for viscous
sintering of ceramic powders. The model describes densification and deformation
during sintering through viscous flow mechanisms, accounting for both
densification-driven and stress-driven deformation.
"""

from pyfem.materials.BaseMaterial import BaseMaterial
from pyfem.materials.MatUtils import hydrostaticStress, transform3To2, transform2To3
from numpy import zeros, ones, dot, exp, array, sqrt
from typing import Tuple


class SkorohodOlevsky(BaseMaterial):
    """
    Skorohod-Olevsky viscous sintering model for ceramics.
    
    This class implements a constitutive model for viscous sintering of ceramic
    powders. The model describes the densification and creep deformation during
    sintering through viscous flow mechanisms.
    
    The model is based on two key concepts:
    
    1. **Densification (volumetric flow)**:
       The rate of densification is driven by sintering stress and hydrostatic
       stress:
       
       dρ/dt = (3ρ/2η_v) * (σ_sint - σ_m)
       
       where:
       - ρ: relative density (current density / theoretical density)
       - η_v: volumetric viscosity (depends on density and temperature)
       - σ_sint: sintering stress (capillary-driven stress from surface tension)
       - σ_m: mean (hydrostatic) stress = (σ₁₁ + σ₂₂ + σ₃₃)/3
    
    2. **Deviatoric deformation (shape change)**:
       Viscous creep under deviatoric stress:
       
       ε̇_dev = s / (2η_s)
       
       where:
       - s: deviatoric stress tensor
       - η_s: shear viscosity
    
    The viscosities are temperature and density dependent:
    - η_v = η₀ * f_v(ρ) * exp(Q/RT)
    - η_s = η₀ * f_s(ρ) * exp(Q/RT)
    
    Required Properties
    -------------------
    eta0 : float
        Reference viscosity at reference temperature (Pa·s).
        Typical values: 1e10 - 1e14 Pa·s for ceramics at sintering temperatures.
    Q : float
        Activation energy for viscous flow (J/mol).
        Typical values: 300000 - 700000 J/mol.
    R : float, optional
        Universal gas constant (J/(mol·K)). Default: 8.314
    T : float
        Temperature (K). For isothermal sintering analysis.
        Typical sintering temperatures: 1400-1800 K for ceramics.
    rho0 : float
        Initial relative density (green density).
        Typical values: 0.5 - 0.7 for powder compacts.
    sigma_sint : float
        Sintering stress (Pa). Related to surface tension and particle size:
        σ_sint ≈ 3γ/r (γ: surface energy, r: particle radius).
        Typical values: 1e5 - 1e7 Pa.
    n_vol : float, optional
        Viscosity exponent for volumetric flow. Default: 2.0
    n_shear : float, optional
        Viscosity exponent for shear flow. Default: 1.0
    
    Examples
    --------
    Material properties for alumina sintering at 1600 K:
    
    <Material>
        type = "SkorohodOlevsky";
        eta0 = 1.0e12;          # Pa·s - reference viscosity
        Q = 500000.0;           # J/mol - activation energy
        T = 1600.0;             # K - sintering temperature
        rho0 = 0.6;             # Initial relative density
        sigma_sint = 1.0e6;     # Pa - sintering stress
        n_vol = 2.0;            # Volumetric viscosity exponent
        n_shear = 1.0;          # Shear viscosity exponent
        R = 8.314;              # J/(mol·K) - gas constant
    </Material>
    
    Notes
    -----
    Physical Interpretation:
    - Relative density ρ increases from initial value (e.g., 0.6) towards 1.0
    - Sintering stress drives densification even without external load
    - Applied compressive stress accelerates densification
    - Applied tensile stress can slow or reverse densification (swelling)
    - Shear deformation occurs through grain boundary sliding
    
    The model captures:
    - Free sintering (densification without applied stress)
    - Pressure-assisted sintering (hot pressing)
    - Creep deformation under stress
    - Temperature-dependent sintering kinetics
    
    References
    ----------
    - Skorohod, V.V. (1972). "Rheological basis of the theory of sintering."
      Naukova Dumka, Kiev.
    - Olevsky, E.A. (1998). "Theory of sintering: from discrete to continuum."
      Materials Science and Engineering: R: Reports, 23(2), 41-100.
    - Olevsky, E.A. and Froyen, L. (2006). "Constitutive modeling of spark-plasma
      sintering of conductive materials." Scripta Materialia, 55(12), 1175-1178.
    """

    def __init__(self, props) -> None:
        """
        Initialize the Skorohod-Olevsky material model.
        
        Parameters
        ----------
        props : Properties
            Material properties object containing eta0, Q, T, rho0, sigma_sint,
            and optional parameters.
        """
        BaseMaterial.__init__(self, props)

        # Set default values
        if not hasattr(self, 'R'):
            self.R = 8.314  # J/(mol·K)
        
        if not hasattr(self, 'n_vol'):
            self.n_vol = 2.0
        
        if not hasattr(self, 'n_shear'):
            self.n_shear = 1.0

        # Validate required parameters
        required_params = ['eta0', 'Q', 'T', 'rho0', 'sigma_sint']
        for param in required_params:
            if not hasattr(self, param):
                raise ValueError(f"Required parameter '{param}' must be specified for SkorohodOlevsky model")

        # Validate parameter ranges
        if self.eta0 <= 0:
            raise ValueError(f"Reference viscosity eta0 must be positive, got {self.eta0}")
        
        if self.Q <= 0:
            raise ValueError(f"Activation energy Q must be positive, got {self.Q}")
        
        if self.T <= 0:
            raise ValueError(f"Temperature T must be positive, got {self.T}")
        
        if not 0 < self.rho0 <= 1.0:
            raise ValueError(f"Initial relative density rho0 must be in (0, 1], got {self.rho0}")

        # Compute temperature-dependent reference viscosity
        # eta_eff = eta0 * exp(Q / RT)
        self.eta_ref = self.eta0 * exp(self.Q / (self.R * self.T))

        # Initialize history variables
        self.setHistoryParameter('rho', self.rho0)  # Current relative density
        self.setHistoryParameter('strain', zeros(6))
        self.setHistoryParameter('strain_visc', zeros(6))  # Viscous strain
        self.setHistoryParameter('sigma', zeros(6))
        self.setHistoryParameter('time_old', 0.0)

        self.commitHistory()

        # Set output labels
        self.outLabels = ["S11", "S22", "S33", "S23", "S13", "S12", "Density", "VolumeStrain"]
        self.outData = zeros(8)

    def getStress(self, kinematics) -> Tuple[array, array]:
        """
        Compute stress and tangent stiffness for the Skorohod-Olevsky model.
        
        The algorithm integrates the viscous sintering equations:
        1. Update relative density from volumetric viscous flow
        2. Compute deviatoric viscous strain increment
        3. Calculate elastic response with density-dependent moduli
        4. Assemble consistent tangent matrix
        
        Parameters
        ----------
        kinematics : Kinematics
            Kinematics object containing strain information.
            Required attributes:
            - dstrain: strain increment (3 or 6 components)
        
        Returns
        -------
        sigma : ndarray
            Stress tensor (3 or 6 components).
        tang : ndarray
            Tangent stiffness matrix (3x3 or 6x6).
        
        Notes
        -----
        The densification rate follows:
        dρ/dt = (3ρ/2η_v) * (σ_sint - σ_m)
        
        where η_v = η_ref * (1/ρ^n_vol - 1)
        
        The deviatoric viscous strain rate follows:
        dε_dev/dt = s / (2η_s)
        
        where η_s = η_ref * (1/ρ^n_shear - 1)
        """
        
        # Get history variables
        strain_old = self.getHistoryParameter('strain')
        strain_visc_old = self.getHistoryParameter('strain_visc')
        rho = self.getHistoryParameter('rho')
        time_old = self.getHistoryParameter('time_old')

        # Get time increment
        time_new = self.solverStat.time
        dtime = time_new - time_old

        # Convert strain to 6 components if needed
        if len(kinematics.dstrain) == 6:
            dstrain = kinematics.dstrain
        else:
            dstrain = transform2To3(kinematics.dstrain)

        # Update total strain
        strain = strain_old + dstrain

        # Compute density-dependent viscosities
        # Skorohod function: η = η_ref * (1/ρ^n - 1)
        if rho < 0.999:  # Avoid singularity at full density
            eta_vol = self.eta_ref * (rho ** (-self.n_vol) - 1.0)
            eta_shear = self.eta_ref * (rho ** (-self.n_shear) - 1.0)
        else:
            # Nearly full density - use very high viscosity (essentially elastic)
            eta_vol = 1.0e20
            eta_shear = 1.0e20

        # Elastic trial stress (based on elastic strain = total - viscous)
        strain_elastic = strain - strain_visc_old
        
        # Compute density-dependent elastic moduli
        # Simple model: E and bulk modulus scale with relative density
        # More sophisticated models use (ρ - ρ0)^n / (1 - ρ0)^n
        rho_factor = rho / self.rho0
        
        # Effective moduli (simplified density dependence)
        # For sintering, typically E ~ ρ^2 to ρ^3
        E_eff = 100.0e9 * (rho_factor ** 2.5)  # Effective Young's modulus
        nu_eff = 0.25  # Assumed constant
        
        # Compute elastic constants with current density
        ebulk3 = E_eff / (1.0 - 2.0 * nu_eff)
        eg2 = E_eff / (1.0 + nu_eff)
        eg = 0.5 * eg2
        elam = (ebulk3 - eg2) / 3.0

        # Construct elastic stiffness matrix
        C = zeros(shape=(6, 6))
        C[:3, :3] = elam
        C[0, 0] += eg2
        C[1, 1] = C[0, 0]
        C[2, 2] = C[0, 0]
        C[3, 3] = eg
        C[4, 4] = C[3, 3]
        C[5, 5] = C[3, 3]

        # Compute trial stress
        sigma_trial = dot(C, strain_elastic)

        # Extract hydrostatic and deviatoric components
        sigma_m = hydrostaticStress(sigma_trial)  # Mean stress
        sigma_dev = sigma_trial.copy()
        sigma_dev[:3] = sigma_dev[:3] - sigma_m * ones(3)  # Deviatoric stress

        # Initialize viscous strain increment
        dstrain_visc = zeros(6)

        if dtime > 0:
            # 1. VOLUMETRIC FLOW (Densification)
            # dρ/dt = (3ρ/2η_v) * (σ_sint - σ_m)
            driving_stress = self.sigma_sint - sigma_m
            drho_dt = (3.0 * rho / (2.0 * eta_vol)) * driving_stress
            
            # Update density with limiter
            drho = drho_dt * dtime
            rho_new = min(rho + drho, 1.0)  # Cannot exceed full density
            rho_new = max(rho_new, self.rho0)  # Cannot decrease below initial
            
            # Volumetric viscous strain (negative for densification)
            # dε_vol = -dρ/ρ (approximate for small changes)
            dstrain_vol_visc = -drho / rho
            dstrain_visc[:3] = dstrain_vol_visc / 3.0  # Distribute to normal components
            
            # Update density
            self.setHistoryParameter('rho', rho_new)

            # 2. DEVIATORIC FLOW (Shape change)
            # dε_dev/dt = s / (2η_s)
            dstrain_dev_visc = sigma_dev / (2.0 * eta_shear) * dtime
            dstrain_visc += dstrain_dev_visc

        # Update viscous strain
        strain_visc = strain_visc_old + dstrain_visc
        self.setHistoryParameter('strain_visc', strain_visc)

        # Recompute elastic strain and stress with updated viscous strain
        strain_elastic = strain - strain_visc
        sigma = dot(C, strain_elastic)

        # Compute tangent stiffness
        # For viscous materials: C_tang = C / (1 + C * Δt / η)
        # Simplified: use elastic stiffness (more complex tangent possible)
        tang = C.copy()
        
        # Add viscous compliance contribution for consistent tangent
        if dtime > 0:
            # Volumetric part
            K = ebulk3 / 3.0  # Bulk modulus
            K_tang = K / (1.0 + K * dtime * 3.0 / (2.0 * eta_vol))
            
            # Shear part
            G = eg  # Shear modulus
            G_tang = G / (1.0 + G * dtime / eta_shear)
            
            # Reassemble tangent with viscous contributions
            lam_tang = K_tang - 2.0 * G_tang / 3.0
            tang[:3, :3] = lam_tang
            tang[0, 0] += 2.0 * G_tang
            tang[1, 1] = tang[0, 0]
            tang[2, 2] = tang[0, 0]
            tang[3, 3] = G_tang
            tang[4, 4] = G_tang
            tang[5, 5] = G_tang

        # Update history variables
        self.setHistoryParameter('sigma', sigma)
        self.setHistoryParameter('strain', strain)
        self.setHistoryParameter('time_old', time_new)

        # Store output data
        self.outData[:6] = sigma
        self.outData[6] = rho
        self.outData[7] = strain_elastic[0] + strain_elastic[1] + strain_elastic[2]  # Volumetric strain

        # Return stress and tangent in appropriate format
        if len(kinematics.dstrain) == 6:
            return sigma, tang
        else:
            return transform3To2(sigma, tang)

    def maximum_principal_stress(self, stress: array) -> float:
        """
        Compute maximum principal stress (alternative to von Mises).
        
        This is more appropriate for brittle materials like ceramics
        that fail primarily under tensile loading.
        
        Parameters
        ----------
        stress : ndarray
            Stress tensor in Voigt notation [σ11, σ22, σ33, σ23, σ13, σ12].
        
        Returns
        -------
        float
            Maximum principal stress (most tensile).
        """
        # For plane stress (σ33 = 0)
        s11, s22 = stress[0], stress[1]
        s12 = stress[5] if len(stress) == 6 else stress[2]
        
        # Principal stresses in 2D
        s_mean = 0.5 * (s11 + s22)
        s_dev = 0.5 * ((s11 - s22) ** 2 + 4 * s12 ** 2) ** 0.5
        
        sigma1 = s_mean + s_dev
        sigma2 = s_mean - s_dev
        
        # Return maximum (most tensile)
        return max(sigma1, sigma2, 0.0)
