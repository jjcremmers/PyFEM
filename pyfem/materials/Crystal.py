# SPDX-License-Identifier: MIT
# Copyright (c) 2011–2026 Joris J.C. Remmers

"""
Crystal plasticity material model for single crystals.

This module implements a crystal plasticity constitutive model based on the
Schmid law for plastic slip in single crystals. The implementation supports:
- Multiple slip systems (up to 3 sets)
- Cubic crystal structures
- Rate-dependent plasticity (power-law viscoplasticity)
- Self- and latent-hardening
- Finite rotation and finite strain (optional)
- Implicit integration with Newton-Raphson iteration

The model is based on the classical crystal plasticity formulation by
Peirce, Shih and Needleman (1984) and incorporates the Bassani and Wu
hardening law with corrections by Kysar (1997).
"""

import numpy as np
from numpy import array, outer, sqrt, abs
from numpy.linalg import norm, solve
from typing import Tuple, List, Any
from pyfem.materials.BaseMaterial import BaseMaterial


class Crystal(BaseMaterial):
    """
    Crystal plasticity constitutive model for single crystals.
    
    This class implements a rate-dependent crystal plasticity model with
    multiple slip systems. The plastic deformation occurs through crystallographic
    slip on specific slip planes in specific slip directions (Schmid's law).
    
    Parameters
    ----------
    props : object
        Properties object containing material parameters:
        
        Elastic properties:
        - E : float, Young's modulus (for isotropic)
        - nu : float, Poisson's ratio (for isotropic)
        - c11, c12, c44 : float, elastic constants (for cubic)
        
        Slip system parameters:
        - nsets : int, number of slip system sets (1-3)
        - plane1, dir1 : array, slip plane normal and direction for set 1
        - plane2, dir2 : array, slip plane normal and direction for set 2 (optional)
        - plane3, dir3 : array, slip plane normal and direction for set 3 (optional)
        
        Crystal orientation:
        - orient1_local, orient1_global : array, first orientation vector
        - orient2_local, orient2_global : array, second orientation vector
        
        Viscoplasticity parameters:
        - gamma0 : float, reference shear strain rate
        - m : float, rate sensitivity exponent
        
        Hardening parameters:
        - tau0 : float, initial critical resolved shear stress
        - h0 : float, initial hardening modulus
        - taus : float, saturation stress
        - a : float, hardening exponent
        - qab : float, latent hardening ratio
        
        Integration parameters:
        - theta : float, implicit integration parameter (0.5 recommended)
        - nlgeom : bool, use finite deformation (default: False)
        - maxiter : int, maximum Newton-Raphson iterations (default: 20)
        - tol : float, convergence tolerance (default: 1e-6)
    """

    def __init__(self, props: Any) -> None:
        """Initialize the crystal plasticity material model."""

        self.theta = 0.5  # Implicit integration parameter
        self.nlgeom = False  # Finite deformation flag
        self.maxiter = 20  # Max iterations
        self.tol = 1e-6  # Convergence tolerance

        BaseMaterial.__init__(self, props)
      
        # Elastic properties
        self._setupElasticProperties()
        
        # Crystal orientation
        self._setupCrystalOrientation()
        
        # Slip systems
        self._setupSlipSystems()
        
        # Viscoplasticity parameters
        self.gamma0 = getattr(props, 'gamma0', 0.001)  # Reference strain rate
        self.m = getattr(props, 'm', 20.0)  # Rate sensitivity exponent
        
        # Hardening parameters
        self.tau0 = getattr(props, 'tau0', 1.0)  # Initial CRSS
        self.h0 = getattr(props, 'h0', 1.0)  # Initial hardening modulus
        self.taus = getattr(props, 'taus', 1.5)  # Saturation stress
        self.a = getattr(props, 'a', 2.0)  # Hardening exponent
        self.qab = getattr(props, 'qab', 1.4)  # Latent hardening ratio
        
        # Initialize history variables
        self._initializeHistory()
        
        # Set output labels
        self._setupOutput()
        
    def _setupElasticProperties(self) -> None:
        """Setup elastic stiffness matrix."""
        
        # Check if cubic crystal constants are provided
        if hasattr(self, 'c11') and hasattr(self, 'c12') and hasattr(self, 'c44'):
            # Cubic crystal
            self.Dlocal = np.zeros((6, 6))
            
            # Diagonal terms (normal stresses)
            self.Dlocal[0:3, 0:3] = self.c12
            np.fill_diagonal(self.Dlocal[0:3, 0:3], self.c11)
            
            # Shear terms
            self.Dlocal[3:6, 3:6] = self.c44 * np.eye(3)
        else:
            # Isotropic material (as fallback)
            E   = self.E
            nu  = self.nu
            lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
            mu  = E / (2.0 * (1.0 + nu))
            
            self.Dlocal = np.zeros((6, 6))
            
            # Diagonal terms (normal stresses)
            self.Dlocal[0:3, 0:3] = lam
            np.fill_diagonal(self.Dlocal[0:3, 0:3], lam + 2.0 * mu)
            
            # Shear terms
            self.Dlocal[3:6, 3:6] = mu * np.eye(3)
    
    def _setupCrystalOrientation(self) -> None:
        """Setup crystal orientation rotation matrix."""
        
        # Get orientation vectors if provided
        if hasattr(self, 'orient1_local') and hasattr(self, 'orient1_global'):
            # Compute rotation matrix from local to global
            v1_local = np.array(self.orient1_local)
            v1_global = np.array(self.orient1_global)
            
            if hasattr(self, 'orient2_local') and hasattr(self, 'orient2_global'):
                v2_local = np.array(self.orient2_local)
                v2_global = np.array(self.orient2_global)
            else:
                # Use default second vector
                v2_local = np.array([0.0, 1.0, 0.0])
                v2_global = np.array([0.0, 1.0, 0.0])
            
            self.R = self._computeRotationMatrix(v1_local, v1_global, 
                                                   v2_local, v2_global)
        else:
            # Identity rotation (no crystal orientation)
            self.R = np.eye(3)
        
        # Transform elastic matrix to global coordinates
        self._transformElasticMatrix()
    
    def _computeRotationMatrix(self, v1_local: np.ndarray, v1_global: np.ndarray, 
                                v2_local: np.ndarray, v2_global: np.ndarray) -> np.ndarray:
        """
        Compute rotation matrix from two non-parallel vectors.
        
        Parameters
        ----------
        v1_local, v1_global : ndarray
            First vector in local and global systems
        v2_local, v2_global : ndarray
            Second vector in local and global systems
            
        Returns
        -------
        R : ndarray (3x3)
            Rotation matrix from local to global coordinates
        """
        
        # Normalize vectors using helper function
        v1_local = self._normalize(v1_local)
        v1_global = self._normalize(v1_global)
        v2_local = self._normalize(v2_local)
        v2_global = self._normalize(v2_global)
        
        # Create local coordinate system
        e1_local = v1_local
        e3_local = np.cross(v1_local, v2_local)
        e3_local = e3_local / norm(e3_local)
        e2_local = np.cross(e3_local, e1_local)
        
        # Create global coordinate system
        e1_global = v1_global
        e3_global = np.cross(v1_global, v2_global)
        e3_global = e3_global / norm(e3_global)
        e2_global = np.cross(e3_global, e1_global)
        
        # Local basis matrix (columns are basis vectors)
        A_local = np.column_stack([e1_local, e2_local, e3_local])
        A_global = np.column_stack([e1_global, e2_global, e3_global])
        
        # Rotation matrix: R = A_global * A_local^T
        R = A_global @ A_local.T
        
        return R
    
    def _transformElasticMatrix(self) -> None:
        """Transform elastic matrix from local to global coordinates."""
        
        # Transformation matrix for 4th order tensor (Voigt notation)
        T = self._getVoigtTransformation(self.R)
        
        # Transform: D_global = T * D_local * T^T
        self.D = T @ self.Dlocal @ T.T
    
    def _getVoigtTransformation(self, R: np.ndarray) -> np.ndarray:
        """
        Get transformation matrix for 4th order tensor in Voigt notation.
        
        Parameters
        ----------
        R : ndarray (3x3)
            Rotation matrix
            
        Returns
        -------
        T : ndarray (6x6)
            Transformation matrix in Voigt notation
        """
        
        T = np.zeros((6, 6))
        
        # Normal-normal components
        T[0:3, 0:3] = R**2
        
        # Shear-shear components using vectorized operations
        shear_indices = [(1, 2), (0, 2), (0, 1)]  # Index pairs for shear components
        for i, (j, k) in enumerate(shear_indices):
            for ii, (jj, kk) in enumerate(shear_indices):
                T[3 + i, 3 + ii] = R[j, jj] * R[k, kk] + R[j, kk] * R[k, jj]
        
        # Normal-shear coupling using vectorized operations
        # T[i, 3:6] components for normal stresses (i=0,1,2)
        T[0:3, 3] = 2.0 * R[0:3, 1] * R[0:3, 2]  # Normal to shear 23
        T[0:3, 4] = 2.0 * R[0:3, 0] * R[0:3, 2]  # Normal to shear 13
        T[0:3, 5] = 2.0 * R[0:3, 0] * R[0:3, 1]  # Normal to shear 12
        
        # T[3:6, j] components for shear stresses (j=0,1,2)
        T[3, 0:3] = R[1, 0:3] * R[2, 0:3]  # Shear 23 to normal
        T[4, 0:3] = R[0, 0:3] * R[2, 0:3]  # Shear 13 to normal
        T[5, 0:3] = R[0, 0:3] * R[1, 0:3]  # Shear 12 to normal
        
        return T
    
    def _setupSlipSystems(self) -> None:
        """Setup slip systems based on crystal structure."""
        
        self.nsets = getattr(self, 'nsets', 1)
        
        self.slipSystems = []
        
        # Setup slip system sets
        for iset in range(self.nsets):
            if iset == 0:
                plane = np.array(getattr(self, 'plane1', [1.0, 1.0, 0.0]))
                direction = np.array(getattr(self, 'dir1', [1.0, 1.0, 1.0]))
            elif iset == 1:
                plane = np.array(getattr(self, 'plane2', [1.0, 0.0, 1.0]))
                direction = np.array(getattr(self, 'dir2', [1.0, 1.0, 1.0]))
            elif iset == 2:
                plane = np.array(getattr(self, 'plane3', [0.0, 1.0, 1.0]))
                direction = np.array(getattr(self, 'dir3', [1.0, 1.0, 1.0]))
            
            # Generate all slip systems for this set
            systems = self._generateSlipSystems(plane, direction)
            self.slipSystems.extend(systems)
        
        self.nslip = len(self.slipSystems)
        
        # Compute Schmid factors (slip deformation tensors)
        self._computeSchmidFactors()
    
    def _generateSlipSystems(self, plane_normal: np.ndarray, 
                              slip_direction: np.ndarray) -> List[Tuple[np.ndarray, np.ndarray]]:
        """
        Generate all slip systems for a given plane and direction.
        
        For cubic crystals, this generates symmetrically equivalent systems.
        
        Parameters
        ----------
        plane_normal : ndarray
            Normal to slip plane
        slip_direction : ndarray
            Slip direction
            
        Returns
        -------
        systems : list
            List of tuples (slip_direction, plane_normal)
        """
        
        systems = []
        
        # Generate symmetrically equivalent systems for cubic crystals
        # This is a simplified version - full implementation would generate
        # all 12 systems for {110}<111>, 24 for {111}<110>, etc.
        
        # Normalize
        n = plane_normal / norm(plane_normal)
        s = slip_direction / norm(slip_direction)
        
        # Transform to global coordinates
        n_global = self.R @ n
        s_global = self.R @ s
        
        # For now, just add the primary system
        # A full implementation would add all symmetrically equivalent systems
        systems.append((s_global, n_global))
        
        # Example: add a few more systems by permutation for {110}<111>
        if abs(abs(n[0]) - abs(n[1])) < 1e-6 and abs(n[2]) < 1e-6:
            # {110} type plane
            perms = [
                ([1, 1, 0], [1, 1, 1]),
                ([1, -1, 0], [1, 1, 1]),
                ([1, 0, 1], [1, 1, 1]),
                ([1, 0, -1], [1, 1, 1]),
                ([0, 1, 1], [1, 1, 1]),
                ([0, 1, -1], [1, 1, 1]),
            ]
            
            for pn, ps in perms[1:]:  # Skip first as it's already added
                n_local = np.array(pn, dtype=float)
                s_local = np.array(ps, dtype=float)
                n_local = self._normalize(n_local)
                s_local = self._normalize(s_local)
                
                n_g = self.R @ n_local
                s_g = self.R @ s_local
                systems.append((s_g, n_g))
        
        return systems
    
    def _computeSchmidFactors(self) -> None:
        """Compute Schmid factors for all slip systems."""
        
        self.schmid = np.zeros((6, self.nslip))
        self.spin = np.zeros((3, self.nslip))
        
        for i, (s, n) in enumerate(self.slipSystems):
            # Schmid tensor: symmetric part of s ⊗ n (Voigt notation)
            # For Voigt notation: [11, 22, 33, 23, 13, 12]
            self.schmid[:, i] = [
                s[0] * n[0],                    # σ11
                s[1] * n[1],                    # σ22
                s[2] * n[2],                    # σ33
                s[1] * n[2] + s[2] * n[1],     # σ23 (factor 1 for strain)
                s[0] * n[2] + s[2] * n[0],     # σ13
                s[0] * n[1] + s[1] * n[0]      # σ12
            ]
            
            # Spin tensor: antisymmetric part of s ⊗ n
            # Components: [ω12, ω31, ω23]
            self.spin[:, i] = 0.5 * np.array([
                s[0] * n[1] - s[1] * n[0],     # ω12
                s[2] * n[0] - s[0] * n[2],     # ω31
                s[1] * n[2] - s[2] * n[1]      # ω23
            ])
    
    def _initializeHistory(self) -> None:
        """Initialize history variables."""
        
        # Current slip system strength
        self.setHistoryParameter('tau_c', self.tau0 * np.ones(self.nslip))
        
        # Shear strain in each slip system
        self.setHistoryParameter('gamma', np.zeros(self.nslip))
        
        # Resolved shear stress in each slip system
        self.setHistoryParameter('tau', np.zeros(self.nslip))
        
        # Cumulative shear strain in each slip system
        self.setHistoryParameter('gamma_cum', np.zeros(self.nslip))
        
        # Total cumulative shear strain
        self.setHistoryParameter('gamma_total', 0.0)
        
        # Stress
        self.setHistoryParameter('stress', np.zeros(6))
        
        # Commit initial state
        self.commitHistory()
    
    def _setupOutput(self) -> None:
        """Setup output labels and data."""
        
        self.outLabels = ['S11', 'S22', 'S33', 'S23', 'S13', 'S12', 
                          'GammaTotal']
        
        # Add slip system outputs
        for i in range(min(self.nslip, 5)):  # Limit output to first 5 systems
            self.outLabels.append(f'Gamma{i+1}')
            self.outLabels.append(f'Tau{i+1}')
        
        self.outData = np.zeros(len(self.outLabels))
    
    def getStress(self, kinematics: Any) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute stress and tangent stiffness matrix.
        
        Parameters
        ----------
        kinematics : object
            Kinematics object containing strain increment (kinematics.dstrain)
            
        Returns
        -------
        stress : ndarray
            Cauchy stress tensor (Voigt notation)
        tang : ndarray (6x6)
            Material tangent stiffness matrix
        """
        
        # Get history variables
        tau_c = self.getHistoryParameter('tau_c')
        gamma = self.getHistoryParameter('gamma')
        tau = self.getHistoryParameter('tau')
        gamma_cum = self.getHistoryParameter('gamma_cum')
        gamma_total = self.getHistoryParameter('gamma_total')
        stress = self.getHistoryParameter('stress')
        
        # Strain increment
        dstrain = kinematics.dstrain
        
        # Time increment (assume unit time if not available)
        dt = self.solverStat.dtime
        
        # Trial elastic stress increment
        dstress_trial = self.D @ dstrain
        stress_trial = stress + dstress_trial
        
        # Compute resolved shear stress on each slip system
        tau_trial = self.schmid.T @ stress_trial
        
        # Shear strain rates (initial guess)
        gamma_dot = self._computeSlipRates(tau_trial, tau_c)
        
        # Implicit integration with Newton-Raphson iteration
        if self.theta > 0.0:
            dgamma, stress, tau_c = self._implicitIntegration(
                stress_trial, tau_trial, tau_c, gamma_dot, dt, dstrain
            )
        else:
            # Explicit integration
            dgamma = gamma_dot * dt
            stress, tau_c = self._explicitIntegration(
                stress_trial, tau_c, dgamma, dt
            )
        
        # Update history variables
        gamma += dgamma
        gamma_cum += np.abs(dgamma)
        gamma_total += np.sum(np.abs(dgamma))
        tau = self.schmid.T @ stress
        
        self.setHistoryParameter('stress', stress)
        self.setHistoryParameter('tau_c', tau_c)
        self.setHistoryParameter('gamma', gamma)
        self.setHistoryParameter('tau', tau)
        self.setHistoryParameter('gamma_cum', gamma_cum)
        self.setHistoryParameter('gamma_total', gamma_total)
        
        # Compute tangent stiffness
        tang = self._computeTangent(stress, tau, tau_c, dgamma, dt)
        
        # Store output data
        self._storeOutput(stress, tau, gamma, gamma_total)
        
        return stress, tang

#-------------------------------------------------------------------
#
#-------------------------------------------------------------------

    def _computeSlipRates(self, tau: np.ndarray, tau_c: np.ndarray) -> np.ndarray:
        """
        Compute shear strain rates using power law.
        
        Parameters
        ----------
        tau : ndarray
            Resolved shear stress in slip systems
        tau_c : ndarray
            Critical resolved shear stress
            
        Returns
        -------
        gamma_dot : ndarray
            Shear strain rates
        """
        
        mask = np.abs(tau_c) > 1e-12
        x = np.where(mask, tau / tau_c, 0.0)
        gamma_dot = self.gamma0 * np.sign(tau) * np.abs(x)**self.m
        
        return gamma_dot
    
    def _computeHardeningMatrix(self, gamma_cum: np.ndarray) -> np.ndarray:
        """
        Compute self- and latent-hardening matrix.
        
        Parameters
        ----------
        gamma_cum : ndarray
            Cumulative shear strain in each slip system
            
        Returns
        -------
        H : ndarray (nslip x nslip)
            Hardening matrix
        """
        
        # Create hardening matrix using broadcasting
        # Start with latent hardening for all entries
        H = self.qab * self.h0 * (1.0 - self.tau0 / self.taus)**self.a * np.ones((self.nslip, self.nslip))
        
        # Overwrite diagonal with self-hardening (q = 1.0)
        h_self = self.h0 * (1.0 - self.tau0 / self.taus)**self.a
        np.fill_diagonal(H, h_self)
        
        return H
    
    def _implicitIntegration(self, stress_trial: np.ndarray, tau_trial: np.ndarray, 
                              tau_c0: np.ndarray, gamma_dot0: np.ndarray, 
                              dt: float, dstrain: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Implicit integration with Newton-Raphson iteration.
        
        Parameters
        ----------
        stress_trial : ndarray
            Trial stress
        tau_trial : ndarray
            Trial resolved shear stress
        tau_c0 : ndarray
            Critical resolved shear stress at start of increment
        gamma_dot0 : ndarray
            Initial shear strain rates
        dt : float
            Time increment
        dstrain : ndarray
            Strain increment
            
        Returns
        -------
        dgamma : ndarray
            Shear strain increments
        stress : ndarray
            Updated stress
        tau_c : ndarray
            Updated critical resolved shear stress
        """
        
        dgamma = gamma_dot0 * dt
        stress = stress_trial.copy()
        tau_c = tau_c0.copy()
        
        for iter in range(self.maxiter):
            # Current resolved shear stress
            tau = self.schmid.T @ stress
            
            # Current shear strain rates
            gamma_dot = self._computeSlipRates(tau, tau_c)
            
            # Residual
            res = dgamma - self.theta * dt * gamma_dot - (1.0 - self.theta) * dt * gamma_dot0
            
            # Check convergence
            if norm(res) < self.tol * self.gamma0 * dt:
                break
            
            # Jacobian
            J = self._computeJacobian(stress, tau, tau_c, dt)
            
            # Newton update
            try:
                ddgamma = solve(J, -res)
            except:
                # If singular, use smaller update
                ddgamma = -0.1 * res
            
            dgamma = dgamma + ddgamma
            
            # Update stress
            dstress_plastic = np.zeros(6)
            for i in range(self.nslip):
                dstress_plastic = dstress_plastic + ddgamma[i] * (self.D @ self.schmid[:, i])
            
            stress = stress - dstress_plastic
            
            # Update hardening
            gamma_cum = self.getHistoryParameter('gamma_cum') + abs(dgamma)
            H = self._computeHardeningMatrix(gamma_cum)
            dtau_c = H @ abs(dgamma)
            tau_c = tau_c0 + dtau_c
        
        return dgamma, stress, tau_c
    
    def _explicitIntegration(self, stress_trial: np.ndarray, tau_c0: np.ndarray, 
                              dgamma: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Explicit integration (forward Euler).
        
        Parameters
        ----------
        stress_trial : ndarray
            Trial stress
        tau_c0 : ndarray
            Critical resolved shear stress at start
        dgamma : ndarray
            Shear strain increments
        dt : float
            Time increment
            
        Returns
        -------
        stress : ndarray
            Updated stress
        tau_c : ndarray
            Updated critical resolved shear stress
        """
        
        # Plastic stress correction
        stress = stress_trial - self.D @ (self.schmid @ dgamma)
        
        # Update hardening
        gamma_cum = self.getHistoryParameter('gamma_cum') + abs(dgamma)
        H = self._computeHardeningMatrix(gamma_cum)
        tau_c = tau_c0 + H @ abs(dgamma)
        
        return stress, tau_c
    
    def _computeJacobian(self, stress: np.ndarray, tau: np.ndarray, 
                          tau_c: np.ndarray, dt: float) -> np.ndarray:
        """
        Compute Jacobian matrix for Newton-Raphson iteration.
        
        Parameters
        ----------
        stress : ndarray
            Current stress
        tau : ndarray
            Resolved shear stress
        tau_c : ndarray
            Critical resolved shear stress
        dt : float
            Time increment
            
        Returns
        -------
        J : ndarray (nslip x nslip)
            Jacobian matrix
        """
        
        J = np.eye(self.nslip)
        
        gamma_cum = self.getHistoryParameter('gamma_cum')
        H = self._computeHardeningMatrix(gamma_cum)
        
        # Compute derivative of slip rate w.r.t. resolved shear stress
        mask = np.abs(tau_c) > 1e-12
        x = np.where(mask, tau / tau_c, 0.0)
        
        mask_nonzero = (np.abs(x) > 1e-12) & mask
        dfdtau = np.where(mask_nonzero, 
                  self.gamma0 * self.m * np.abs(x)**(self.m - 1.0) / tau_c,
                  0.0)
        
        # Compute stress change contributions: dtau_dgamma[i,j] = -schmid[i]^T D schmid[j]
        dtau_dgamma = -(self.schmid.T @ self.D @ self.schmid)
        
        # Compute hardening contributions: dtauc_dgamma[i,j] = H[i,j] * sign(tau[j])
        dtauc_dgamma = H * np.sign(tau)[np.newaxis, :]
        
        # Compute ratio tau/tau_c safely
        tau_ratio = np.where(mask, tau / tau_c, 0.0)
        
        # Compute full Jacobian using broadcasting
        # J[i,j] -= theta * dt * dfdtau[i] * (dtau_dgamma[i,j] - tau_ratio[i] * dtauc_dgamma[i,j])
        J -= self.theta * dt * dfdtau[:, np.newaxis] * (dtau_dgamma - tau_ratio[:, np.newaxis] * dtauc_dgamma)
        
        return J
    
    def _computeTangent(self, stress: np.ndarray, tau: np.ndarray, 
                         tau_c: np.ndarray, dgamma: np.ndarray, dt: float) -> np.ndarray:
        """
        Compute consistent tangent stiffness matrix.
        
        Parameters
        ----------
        stress : ndarray
            Current stress
        tau : ndarray
            Resolved shear stress
        tau_c : ndarray
            Critical resolved shear stress
        dgamma : ndarray
            Shear strain increments
        dt : float
            Time increment
            
        Returns
        -------
        tang : ndarray (6x6)
            Tangent stiffness matrix
        """
        
        # Check if plastic
        is_plastic = any(abs(dgamma) > self.tol * self.gamma0 * dt)
        
        if not is_plastic:
            # Elastic tangent
            return self.D
        
        # Compute consistent tangent (algorithmic tangent)
        J = self._computeJacobian(stress, tau, tau_c, dt)
        
        try:
            Jinv = solve(J, np.eye(self.nslip))
        except:
            # Use elastic tangent if Jacobian is singular
            return self.D
        
        # Consistent tangent
        tang = self.D.copy()
        
        # Compute D @ schmid for all slip systems at once
        DS = self.D @ self.schmid  # (6, nslip)
        
        # Compute consistent tangent using broadcasting
        # tang -= sum_ij Jinv[i,j] * (D @ schmid_i) @ (schmid_j^T @ D)
        # = DS @ Jinv @ schmid^T @ D
        tang = tang - DS @ Jinv @ self.schmid.T @ self.D
        
        return tang
    
    def _storeOutput(self, stress: np.ndarray, tau: np.ndarray, 
                      gamma: np.ndarray, gamma_total: float) -> None:
        """Store output data for visualization."""
        
        self.outData[0:6] = stress
        self.outData[6] = gamma_total
        
        idx = 7
        for i in range(min(self.nslip, 5)):
            self.outData[idx] = gamma[i]
            self.outData[idx + 1] = tau[i]
            idx += 2

    def _normalize(self, v: np.ndarray) -> np.ndarray:
        """
        Normalize a vector.
        
        Parameters
        ----------
        v : ndarray
            Input vector
            
        Returns
        -------
        v_norm : ndarray
            Normalized vector (unit length)
        """
        n = norm(v)
        if n < 1e-12:
            raise ValueError("Cannot normalize zero vector")
        return v / n        
