# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""
Rate-dependent viscoplastic material model using Perzyna formulation.

This module implements a viscoplastic constitutive model based on the Perzyna
overstress theory. The material exhibits rate-dependent plastic flow when the
stress exceeds the yield surface, with the viscoplastic strain rate proportional
to the overstress.
"""

from pyfem.materials.BaseMaterial import BaseMaterial
from pyfem.materials.MatUtils import vonMisesStress, hydrostaticStress, transform3To2, transform2To3
from numpy import zeros, ones, dot, outer, array, sqrt
from typing import Tuple


class ViscoPlasticity(BaseMaterial):
    """
    Rate-dependent viscoplastic material model using Perzyna formulation.
    
    This class implements a viscoplastic constitutive model where plastic flow
    is rate-dependent. The model uses the Perzyna overstress formulation where
    the viscoplastic strain rate is given by:
    
    d(εᵖ)/dt = γ * <Φ(f)>ⁿ * ∂f/∂σ
    
    where:
    - γ: fluidity parameter (inverse of viscosity)
    - Φ(f): overstress function = (f - f_yield) / f_yield
    - n: rate sensitivity exponent
    - f: von Mises equivalent stress
    - <>: Macaulay brackets (returns value if positive, zero otherwise)
    
    The model combines:
    1. Elastic behavior below yield stress
    2. Rate-dependent plastic flow above yield stress
    3. Optional strain hardening
    
    Required Properties
    -------------------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.
    syield : float
        Initial yield stress (von Mises).
    gamma : float
        Fluidity parameter (1/viscosity), controls rate sensitivity.
        Typical range: 1e-6 to 1e-2 for metals, higher for polymers.
    n : float, optional
        Rate sensitivity exponent. Default: 1.0 (linear viscosity).
        n > 1 increases rate sensitivity at high overstress.
    hard : float, optional
        Linear hardening modulus. Default: 0.0 (perfectly plastic).
    
    Examples
    --------
    Material properties in input file:
    
    <Material>
        type = "ViscoPlasticity";
        E = 200000.0;
        nu = 0.3;
        syield = 250.0;
        gamma = 0.001;
        n = 1.0;
        hard = 1000.0;
    </Material>
    
    Notes
    -----
    - For γ → 0, the model approaches rate-independent plasticity
    - For γ → ∞, the model approaches perfect viscous behavior
    - The time integration uses a backward Euler scheme for stability
    """

    def __init__(self, props) -> None:
        """
        Initialize the viscoplastic material model.
        
        Parameters
        ----------
        props : Properties
            Material properties object containing E, nu, syield, gamma, and
            optional parameters n and hard.
        """
        BaseMaterial.__init__(self, props)

        print(self)
        # Set default values
        if not hasattr(self, 'n'):
            self.n = 1.0
        
        if not hasattr(self, 'hard'):
            self.hard = 0.0

        # Validate required parameters
        if not hasattr(self, 'gamma'):
            raise ValueError("Fluidity parameter 'gamma' must be specified for ViscoPlasticity")
        
        if not hasattr(self, 'syield'):
            raise ValueError("Yield stress 'syield' must be specified for ViscoPlasticity")

        # Compute elastic constants
        self.ebulk3 = self.E / (1.0 - 2.0 * self.nu)
        self.eg2 = self.E / (1.0 + self.nu)
        self.eg = 0.5 * self.eg2
        self.eg3 = 3.0 * self.eg
        self.elam = (self.ebulk3 - self.eg2) / 3.0

        # Construct elastic stiffness matrix
        self.ctang = zeros(shape=(6, 6))
        self.ctang[:3, :3] = self.elam

        self.ctang[0, 0] += self.eg2
        self.ctang[1, 1] = self.ctang[0, 0]
        self.ctang[2, 2] = self.ctang[0, 0]

        self.ctang[3, 3] = self.eg
        self.ctang[4, 4] = self.ctang[3, 3]
        self.ctang[5, 5] = self.ctang[3, 3]

        # Initialize history variables
        self.setHistoryParameter('sigma', zeros(6))
        self.setHistoryParameter('eelas', zeros(6))
        self.setHistoryParameter('eplas', zeros(6))
        self.setHistoryParameter('eqplas', 0.0)
        self.setHistoryParameter('time_old', 0.0)

        self.commitHistory()

        # Set output labels
        self.outLabels = ["S11", "S22", "S33", "S23", "S13", "S12", "Epl", "EqPl"]
        self.outData = zeros(8)

        # Convergence tolerance
        self.tolerance = 1.0e-8
        self.maxIter = 20

    def getStress(self, kinematics) -> Tuple[array, array]:
        """
        Compute stress and tangent stiffness for the viscoplastic model.
        
        Uses backward Euler time integration to solve the rate-dependent plastic
        flow equations. The algorithm includes:
        1. Elastic predictor step
        2. Check for plastic admissibility
        3. Viscoplastic corrector with local Newton iteration
        
        Parameters
        ----------
        kinematics : Kinematics
            Kinematics object containing strain increment.
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
        The viscoplastic strain increment is computed using:
        Δεᵖ = Δt * γ * <(σ_trial - σ_yield)/σ_yield>ⁿ * N
        
        where N is the flow direction (normal to yield surface).
        """
        
        # Get history variables
        eelas = self.getHistoryParameter('eelas')
        eplas = self.getHistoryParameter('eplas')
        eqplas = self.getHistoryParameter('eqplas')
        sigma = self.getHistoryParameter('sigma')
        time_old = self.getHistoryParameter('time_old')

        # Get time increment
        time_new = self.solverStat.time
        dtime = time_new - time_old

        # Convert strain to 6 components if needed
        if len(kinematics.dstrain) == 6:
            dstrain = kinematics.dstrain
        else:
            dstrain = transform2To3(kinematics.dstrain)

        # Elastic predictor
        eelas_trial = eelas + dstrain
        sigma_trial = dot(self.ctang, eelas_trial)

        # Compute von Mises stress
        smises = vonMisesStress(sigma_trial)

        # Current yield stress including hardening
        syield_current = self.syield + self.hard * eqplas

        # Initialize tangent with elastic stiffness
        tang = self.ctang.copy()

        # Check for viscoplastic loading
        if smises > syield_current and dtime > 0:
            # Compute overstress ratio
            overstress = (smises - syield_current) / syield_current

            # Viscoplastic multiplier (per unit time)
            gamma_eff = self.gamma * (overstress ** self.n)

            # Incremental viscoplastic strain
            deqpl = gamma_eff * dtime

            # Extract deviatoric stress and compute flow direction
            shydro = hydrostaticStress(sigma_trial)
            flow = sigma_trial.copy()
            flow[:3] = flow[:3] - shydro * ones(3)
            flow *= 1.0 / smises

            # Local Newton iteration for viscoplastic corrector
            converged = False
            for iter in range(self.maxIter):
                # Current yield stress
                syield_iter = self.syield + self.hard * (eqplas + deqpl)

                # Residual
                residual = smises - self.eg3 * deqpl - syield_iter

                # Check convergence
                if abs(residual) < self.tolerance * self.syield:
                    converged = True
                    break

                # Jacobian (derivative of residual w.r.t. deqpl)
                jacobian = -self.eg3 - self.hard

                # Newton update
                deqpl_inc = -residual / jacobian
                deqpl += deqpl_inc

            if not converged:
                import warnings
                warnings.warn(f"Viscoplastic iteration did not converge after {self.maxIter} iterations")

            # Update plastic strain
            eplas[:3] += 1.5 * flow[:3] * deqpl
            eplas[3:] += 3.0 * flow[3:] * deqpl

            # Update elastic strain
            eelas[:3] = eelas_trial[:3] - 1.5 * flow[:3] * deqpl
            eelas[3:] = eelas_trial[3:] - 3.0 * flow[3:] * deqpl

            # Update equivalent plastic strain
            eqplas += deqpl

            # Compute final stress
            syield_final = self.syield + self.hard * eqplas
            sigma = flow * syield_final
            sigma[:3] += shydro * ones(3)

            # Compute consistent tangent (algorithmic tangent)
            # Effective shear moduli
            effg = self.eg * syield_final / smises
            effg2 = 2.0 * effg
            effg3 = 3.0 * effg
            efflam = (self.ebulk3 - effg2) / 3.0

            # Rate-dependent hardening contribution
            rate_factor = self.gamma * self.n * (overstress ** (self.n - 1)) * dtime / syield_current
            effhdr = self.eg3 * (self.hard + rate_factor) / (self.eg3 + self.hard + rate_factor) - effg3

            # Construct tangent
            tang = zeros(shape=(6, 6))
            tang[:3, :3] = efflam

            for i in range(3):
                tang[i, i] += effg2
                tang[i + 3, i + 3] += effg

            # Add hardening contribution
            tang += effhdr * outer(flow, flow)

        else:
            # Elastic response
            eelas = eelas_trial
            sigma = sigma_trial

        # Update history variables
        self.setHistoryParameter('eelas', eelas)
        self.setHistoryParameter('eplas', eplas)
        self.setHistoryParameter('sigma', sigma)
        self.setHistoryParameter('eqplas', eqplas)
        self.setHistoryParameter('time_old', time_new)

        # Store output data
        self.outData[:6] = sigma
        self.outData[6] = eplas[0]
        self.outData[7] = eqplas

        # Return stress and tangent in appropriate format
        if len(kinematics.dstrain) == 6:
            return sigma, tang
        else:
            return transform3To2(sigma, tang)
