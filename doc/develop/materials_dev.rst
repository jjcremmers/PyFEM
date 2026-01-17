=============================
Developing Material Models
=============================

This guide explains how to implement new material models in PyFEM. Material
models (constitutive laws) define the relationship between stress and strain,
forming the core of finite element analysis for solid mechanics problems.

Overview
--------

Material models in PyFEM are responsible for:

- Computing stress tensors from strain tensors
- Providing tangent moduli (material stiffness) for Newton-Raphson iteration
- Managing internal variables for history-dependent materials (plasticity, damage)
- Providing output quantities for post-processing

All material models inherit from the ``BaseMaterial`` class and implement
methods that are called by element formulations during assembly.

Material Class Structure
------------------------

Base Class
~~~~~~~~~~

All material models must inherit from ``BaseMaterial`` located in
``pyfem/materials/BaseMaterial.py``. The base class handles:

- Property management (E, nu, etc.)
- Output data storage
- Common utility functions

Required Methods
~~~~~~~~~~~~~~~~

A material model must implement:

.. code-block:: python

   class MyMaterial(BaseMaterial):
   
       def __init__(self, props):
           """Initialize material with properties."""
           BaseMaterial.__init__(self, props)
           # Initialize material parameters and state
           
       def getStress(self, deformation):
           """Compute stress and tangent modulus.
           
           Args:
               deformation: Kinematics object with strain data
               
           Returns:
               tuple: (sigma, tangent) - stress vector and tangent matrix
           """
           # Implementation here
           return sigma, tangent

Optional Methods
~~~~~~~~~~~~~~~~

For advanced features:

.. code-block:: python

   def reset(self):
       """Reset internal variables for new load step."""
       pass
       
   def commit(self):
       """Commit current state after convergence."""
       pass

Implementation Examples
-----------------------

Example 1: Elastic Material (Plane Stress)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example implements a plane stress elastic material following Chapter 3
of the book *"Non-Linear Finite Element Analysis of Solids and Structures"*
by de Borst et al.

The plane stress constitutive relation (equation 3.26) is:

.. math::

   \begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{12} \end{bmatrix}
   = \frac{E}{1-\nu^2}
   \begin{bmatrix} 
   1 & \nu & 0 \\
   \nu & 1 & 0 \\
   0 & 0 & \frac{1-\nu}{2}
   \end{bmatrix}
   \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{22} \\ \gamma_{12} \end{bmatrix}

.. code-block:: python

   # SPDX-License-Identifier: MIT
   # Copyright (c) 2011–2026 Your Name

   from pyfem.materials.BaseMaterial import BaseMaterial
   import numpy as np

   class PlaneStress(BaseMaterial):
       """Plane stress elastic material model.
       
       Implements the plane stress assumption where σ_33 = 0.
       Based on equation (3.26) in de Borst et al.
       
       Properties:
           E: Young's modulus
           nu: Poisson's ratio
       """
   
       def __init__(self, props):
           """Initialize plane stress material.
           
           Args:
               props: Properties object containing E and nu
           """
           BaseMaterial.__init__(self, props)
           
           # Build elasticity matrix H (equation 3.26)
           self.H = np.zeros((3, 3))
           
           factor = self.E / (1.0 - self.nu * self.nu)
           
           self.H[0, 0] = factor
           self.H[0, 1] = factor * self.nu
           self.H[1, 0] = factor * self.nu
           self.H[1, 1] = factor
           self.H[2, 2] = self.E / (2.0 * (1.0 + self.nu))
           
           # Define output labels
           self.outLabels = ["S11", "S22", "S12"]
   
       def getStress(self, deformation):
           """Compute stress from strain.
           
           Args:
               deformation: Kinematics object with strain vector
                           [ε_11, ε_22, γ_12]
           
           Returns:
               tuple: (sigma, H) where
                   sigma: Stress vector [σ_11, σ_22, σ_12]
                   H: Tangent elasticity matrix (constant)
           """
           # Linear elastic: σ = H ε
           sigma = self.H @ deformation.strain
           
           # Store output data for post-processing
           self.outData = sigma
           
           # Return stress and tangent (same as H for elastic)
           return sigma, self.H
   
       def getTangent(self):
           """Return tangent modulus matrix.
           
           For elastic materials, this is constant.
           """
           return self.H

Example 2: Von Mises Plasticity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example implements J2 plasticity with isotropic hardening, following
Chapter 6 of the book. The implementation uses a return mapping algorithm
(Box 6.1).

.. code-block:: python

   from pyfem.materials.BaseMaterial import BaseMaterial
   import numpy as np
   from numpy.linalg import norm

   class VonMises(BaseMaterial):
       """Von Mises plasticity with isotropic hardening.
       
       Implements J2 plasticity following Chapter 6 of de Borst et al.
       Uses return mapping algorithm (Box 6.1, page 185).
       
       Properties:
           E: Young's modulus
           nu: Poisson's ratio
           sY: Initial yield stress
           hard: Hardening modulus
       """
   
       def __init__(self, props):
           """Initialize von Mises material."""
           BaseMaterial.__init__(self, props)
           
           # Elastic properties
           K = self.E / (3.0 * (1.0 - 2.0 * self.nu))  # Bulk modulus
           G = self.E / (2.0 * (1.0 + self.nu))        # Shear modulus
           
           self.K = K
           self.G = G
           
           # Plasticity properties
           self.sY = self.sY      # Yield stress (from props)
           self.hard = self.hard  # Hardening modulus (from props)
           
           # Internal variables
           self.epse = 0.0  # Equivalent plastic strain
           self.epse_old = 0.0
           
           self.outLabels = ["S11", "S22", "S33", "S23", "S13", "S12",
                           "eqps"]  # Equivalent plastic strain
   
       def getStress(self, deformation):
           """Compute stress using return mapping algorithm.
           
           Args:
               deformation: Kinematics with strain vector
           
           Returns:
               tuple: (sigma, tangent)
           """
           strain = deformation.strain
           
           # Split strain into volumetric and deviatoric parts
           eps_v = strain[0] + strain[1] + strain[2]  # Volumetric strain
           eps_dev = strain.copy()
           eps_dev[0] -= eps_v / 3.0
           eps_dev[1] -= eps_v / 3.0
           eps_dev[2] -= eps_v / 3.0
           
           # Elastic predictor (equation 6.10)
           p_trial = self.K * eps_v  # Pressure
           s_trial = 2.0 * self.G * eps_dev  # Deviatoric stress
           
           # Compute von Mises equivalent stress (equation 6.6)
           q_trial = np.sqrt(1.5 * (s_trial[0]**2 + s_trial[1]**2 + 
                                     s_trial[2]**2 + 2.0*s_trial[3]**2 +
                                     2.0*s_trial[4]**2 + 2.0*s_trial[5]**2))
           
           # Check yield condition (equation 6.7)
           f_trial = q_trial - (self.sY + self.hard * self.epse_old)
           
           if f_trial <= 0:
               # Elastic step
               sigma = np.zeros(6)
               sigma[0] = s_trial[0] + p_trial
               sigma[1] = s_trial[1] + p_trial
               sigma[2] = s_trial[2] + p_trial
               sigma[3:6] = s_trial[3:6]
               
               self.epse = self.epse_old
               tangent = self.getElasticTangent()
               
           else:
               # Plastic step - return mapping (Box 6.1)
               Dgamma = f_trial / (3.0 * self.G + self.hard)
               
               # Update equivalent plastic strain
               self.epse = self.epse_old + Dgamma
               
               # Return mapping: project to yield surface
               factor = 1.0 - (3.0 * self.G * Dgamma) / q_trial
               s = factor * s_trial
               
               # Reconstruct stress tensor
               sigma = np.zeros(6)
               sigma[0] = s[0] + p_trial
               sigma[1] = s[1] + p_trial
               sigma[2] = s[2] + p_trial
               sigma[3:6] = s[3:6]
               
               # Compute elastoplastic tangent (equation 6.32)
               tangent = self.getPlasticTangent(s, q_trial, Dgamma)
           
           # Store output
           self.outData = np.append(sigma, self.epse)
           
           return sigma, tangent
   
       def getElasticTangent(self):
           """Build elastic tangent matrix."""
           D = np.zeros((6, 6))
           
           # Bulk and shear contributions
           lam = self.K - 2.0 * self.G / 3.0
           
           D[0:3, 0:3] = lam
           D[0, 0] = D[1, 1] = D[2, 2] = lam + 2.0 * self.G
           D[3, 3] = D[4, 4] = D[5, 5] = self.G
           
           return D
   
       def getPlasticTangent(self, s, q, Dgamma):
           """Build consistent elastoplastic tangent (equation 6.32)."""
           D_ep = self.getElasticTangent()
           
           if q > 1e-10:
               # Compute direction tensor n
               n = s / q
               
               # Consistency parameter
               theta = 1.0 - (3.0 * self.G * Dgamma) / q
               
               # Modification for plastic loading
               factor = (6.0 * self.G**2) / (3.0 * self.G + self.hard)
               factor *= (1.0 / q - theta / (q**2))
               
               # Rank-one update (equation 6.32)
               for i in range(6):
                   for j in range(6):
                       D_ep[i, j] -= factor * s[i] * s[j]
           
           return D_ep
   
       def commit(self):
           """Commit converged state."""
           self.epse_old = self.epse
   
       def reset(self):
           """Reset to last converged state."""
           self.epse = self.epse_old

Example 3: Cohesive Zone Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Interface elements require traction-separation laws. Here's an example
implementing the Xu-Needleman model (Chapter 10):

.. code-block:: python

   from pyfem.materials.BaseMaterial import BaseMaterial
   import numpy as np

   class XuNeedleman(BaseMaterial):
       """Xu-Needleman cohesive zone model.
       
       Implements exponential traction-separation law for mode I fracture.
       Based on Chapter 10 of de Borst et al.
       
       Properties:
           Tult: Maximum traction
           Gc: Fracture energy (area under T-δ curve)
       """
   
       def __init__(self, props):
           """Initialize cohesive material."""
           BaseMaterial.__init__(self, props)
           
           # Characteristic opening displacement
           self.delta_c = np.e * self.Gc / self.Tult
           
           self.outLabels = ["Tn", "Tt", "un", "ut"]
   
       def getStress(self, deformation):
           """Compute traction from displacement jump.
           
           Args:
               deformation: Contains displacement jump [Δn, Δt]
           
           Returns:
               tuple: (traction, tangent)
           """
           # Displacement jump
           dn = deformation.strain[0]  # Normal opening
           dt = deformation.strain[1]  # Tangential slip
           
           # Effective opening (mixed mode)
           delta_eff = np.sqrt(dn**2 + dt**2)
           
           if delta_eff < 1e-10:
               # No opening: elastic interface
               return np.zeros(2), self.getInitialStiffness()
           
           # Traction law: T = T_ult * (δ/δ_c) * exp(1 - δ/δ_c)
           ratio = delta_eff / self.delta_c
           T_eff = self.Tult * ratio * np.exp(1.0 - ratio)
           
           # Decompose into normal and tangential
           Tn = T_eff * dn / delta_eff
           Tt = T_eff * dt / delta_eff
           
           traction = np.array([Tn, Tt])
           
           # Compute tangent (derivative of T with respect to δ)
           tangent = self.computeTangent(dn, dt, delta_eff, ratio)
           
           # Store output
           self.outData = np.array([Tn, Tt, dn, dt])
           
           return traction, tangent
   
       def computeTangent(self, dn, dt, delta_eff, ratio):
           """Compute tangent stiffness for cohesive interface."""
           if delta_eff < 1e-10:
               return self.getInitialStiffness()
           
           # Derivative of effective traction
           dT_ddelta = (self.Tult / self.delta_c) * (1.0 - ratio) * \
                       np.exp(1.0 - ratio)
           
           # Build tangent matrix
           K = np.zeros((2, 2))
           
           K[0, 0] = dT_ddelta * (dn / delta_eff)**2
           K[0, 1] = dT_ddelta * dn * dt / (delta_eff**2)
           K[1, 0] = K[0, 1]
           K[1, 1] = dT_ddelta * (dt / delta_eff)**2
           
           return K
   
       def getInitialStiffness(self):
           """Penalty stiffness for small openings."""
           K_penalty = self.Tult / (0.01 * self.delta_c)
           return K_penalty * np.eye(2)

State Models vs. Constitutive Models
-------------------------------------

PyFEM uses a hierarchical material system:

State Models
~~~~~~~~~~~~

State models (like PlaneStress, PlaneStrain) adapt 3D constitutive laws to
specific stress states. They wrap underlying constitutive models.

.. code-block:: python

   class PlaneStrain(BaseMaterial):
       """Wrapper for plane strain condition."""
       
       def __init__(self, props):
           BaseMaterial.__init__(self, props)
           
           # Create underlying 3D model if specified
           if hasattr(props, 'model'):
               self.model = self.create_material(props.model)
           else:
               # Use linear elastic by default
               self.model = None

Nested Usage
~~~~~~~~~~~~

In input files, users can nest constitutive models:

.. code-block:: text

   material = 
   {
     type = "PlaneStrain";
     E    = 210.0e3;
     nu   = 0.3;
     
     model = 
     {
       type = "VonMises";
       sY   = 250.0;
       hard = 1000.0;
     };
   };

Registration and Usage
----------------------

File Location
~~~~~~~~~~~~~

Place your material class in:

- ``pyfem/materials/MyMaterial.py``

Import in __init__.py
~~~~~~~~~~~~~~~~~~~~~

Add to ``pyfem/materials/__init__.py``:

.. code-block:: python

   from .MyMaterial import MyMaterial
   
   __all__ = [
       'MyMaterial',
       # ... other materials
   ]

Using in Input Files
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   ElementGroup = 
   {
     type = "SmallStrainContinuum";
     
     material = 
     {
       type = "MyMaterial";
       E    = 210.0e3;
       nu   = 0.3;
       # ... custom properties
     };
   };

Testing and Validation
----------------------

Unit Tests
~~~~~~~~~~

Test individual material responses:

.. code-block:: python

   import unittest
   from pyfem.materials.MyMaterial import MyMaterial
   
   class TestMyMaterial(unittest.TestCase):
   
       def test_uniaxial_tension(self):
           """Test uniaxial stress state."""
           # Create material
           # Apply strain
           # Check stress matches theory
           pass

Analytical Benchmarks
~~~~~~~~~~~~~~~~~~~~~

Compare against closed-form solutions:

1. **Uniaxial tension/compression**
2. **Pure shear**
3. **Hydrostatic compression**
4. **Cyclic loading** (for plasticity)

Patch Tests
~~~~~~~~~~~

Integrate with element tests to verify consistency.

Best Practices
--------------

1. **Symmetric tangent**: Ensure tangent matrix is symmetric
2. **Positive definiteness**: Check for negative eigenvalues
3. **Consistent units**: Document expected units
4. **Handle singularities**: Check for zero denominators
5. **State management**: Properly commit/reset internal variables
6. **Output data**: Provide useful post-processing quantities
7. **Documentation**: Reference equations from the book

Common Pitfalls
---------------

1. **Forgetting to update internal variables** in commit()
2. **Incorrect tangent modulus** leading to convergence issues
3. **Sign conventions** for stress/strain
4. **Voigt notation** ordering [σ11, σ22, σ33, σ23, σ13, σ12]
5. **Engineering vs. tensorial shear strain** (factor of 2)

References
----------

The theoretical foundation for material models can be found in:

*"Non-Linear Finite Element Analysis of Solids and Structures"*
by R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
John Wiley & Sons, 2012, ISBN 978-0470666449

Key chapters:

- Chapter 3: Constitutive Models
- Chapter 6: Plasticity
- Chapter 7: Damage Mechanics
- Chapter 10: Discontinuities and Localization

See Also
--------

- :doc:`elements_dev` - Implementing element formulations
- :doc:`solvers_dev` - Implementing solution algorithms
- :doc:`io_dev` - Implementing I/O modules
- :doc:`../materials/overview` - Available material models
