========================
FiniteStrainContinuum
========================

The ``FiniteStrainContinuum`` element is a displacement-based finite element 
for 2D and 3D large deformation analysis in continuum mechanics. This element 
formulation accounts for geometric nonlinearity through finite strain kinematics, 
making it suitable for problems involving large displacements, large rotations, 
and moderate to large strains.

The element supports both Total Lagrangian (TL) and Updated Lagrangian (UL) 
formulations, automatically adapts to 2D or 3D analysis based on the problem 
dimension (``rank``), and can be combined with various material models including 
hyperelastic and inelastic constitutive laws.

--------
Overview
--------

Element type: ``FiniteStrainContinuum``

The element implements:

- **Total Lagrangian (TL) formulation**: Reference configuration-based approach 
  using the initial geometry throughout the analysis
- **Updated Lagrangian (UL) formulation**: Current configuration-based approach 
  where the reference is updated at each increment
- **2D formulation**: Plane stress or plane strain analysis with triangular 
  or quadrilateral elements (displacement DOFs: ``u``, ``v``)
- **3D formulation**: Full 3D analysis with tetrahedral or hexahedral elements 
  (displacement DOFs: ``u``, ``v``, ``w``)
- **Geometric nonlinearity**: Full consideration of large displacement and 
  rotation effects through geometric stiffness contributions
- **Material nonlinearity**: Compatible with hyperelastic, plastic, and other 
  nonlinear constitutive models

The element computes the tangent stiffness matrix including both material and 
geometric contributions, internal force vector, and mass matrix for dynamic 
analysis. The formulation includes strain-displacement relations appropriate 
for finite deformations.

----------
Parameters
----------

Mandatory parameters
~~~~~~~~~~~~~~~~~~~~

``type``
  Must be set to ``"FiniteStrainContinuum"``

``material``
  A material block defining the constitutive behavior for finite strains. 
  The material block must include a ``type`` parameter specifying the material 
  model.

  Common material types include:

  - ``"PlaneStress"``: For 2D plane stress analysis (can handle finite strains)
  - ``"PlaneStrain"``: For 2D plane strain analysis
  - ``"Isotropic"``: For 3D isotropic elasticity
  - ``"NeoHookean"``: Hyperelastic material model for rubber-like behavior
  - Other finite strain material models as available

  Within the material block, material-specific parameters must be provided:

  - ``E``: Young's modulus
  - ``nu``: Poisson's ratio
  - Additional parameters depending on the material type

Optional parameters
~~~~~~~~~~~~~~~~~~~

``method``
  Specifies the Lagrangian formulation to use. Options are:

  - ``"TL"``: Total Lagrangian formulation (default)
  - ``"UL"``: Updated Lagrangian formulation

  If not specified, the Total Lagrangian formulation is used.

--------
Examples
--------

Example 1: 2D Nonlinear Beam (Total Lagrangian)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A cantilever beam subjected to large displacements using the default Total 
Lagrangian formulation::

  ContElem = 
  {
    type = "FiniteStrainContinuum";

    material =
    {
      type = "PlaneStress";
      E    = 1.e6;
      nu   = 0.25;
    };
  };

This configuration is used with a nonlinear solver in: 
``examples/ch03/NewtonRaphson.pro``

Example 2: 3D Finite Strain Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 3D structure undergoing large deformations::

  ContElem = 
  {
    type = "FiniteStrainContinuum";

    material =
    {
      type = "Isotropic";
      E    = 1.e6;
      nu   = 0.25;
    };
  };

This configuration is used in: ``examples/ch03/NewtonRaphson3D.pro``

Example 3: Dynamic Analysis with Finite Strains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The element can be used for explicit dynamic analysis with finite strain effects::

  ContElem =
  {
    type = "FiniteStrainContinuum";

    material =
    {
      type = "PlaneStress";
      E    = 1.e6;
      nu   = 0.25;
      rho  = 1000.;
    };
  };

This configuration is used in dynamic stress wave propagation: 
``examples/ch05/StressWave_20x20.pro``

------------------
Additional Examples
------------------

The ``FiniteStrainContinuum`` element is used in various examples demonstrating 
large deformation behavior:

- **Nonlinear beam analysis**: ``examples/ch03/cantilever8.pro``, 
  ``examples/ch03/cantilever8PrescribedDisp.pro``
- **Newton-Raphson examples**: ``examples/ch03/NewtonRaphson.pro``, 
  ``examples/ch03/NewtonRaphson3D.pro``
- **Dynamic analysis**: ``examples/ch05/StressWave_20x20.pro``
- **Contact mechanics**: ``examples/contact/contact_test.pro``, 
  ``examples/contact/contact_test02.pro``
- **Buckling analysis**: ``examples/solver/dissipatedEnergySolver/delam_buckling100.pro``, 
  ``examples/solver/dissipatedEnergySolver/delam_buckling200.pro``

---------------
Solver Requirements
---------------

The ``FiniteStrainContinuum`` element requires a nonlinear solver due to the 
geometric nonlinearity. Use solvers such as:

- ``NonlinearSolver``: Standard Newton-Raphson solution procedure
- ``RiksSolver``: Arc-length method for tracing equilibrium paths
- ``DissipatedEnergySolver``: For problems with softening behavior
- ``ExplicitSolver``: For explicit dynamic analysis

See Also
--------

- :doc:`smallstraincontinuum` - Small strain continuum element
- :doc:`materials` - Available material models
- :doc:`solvers` - Nonlinear solution procedures
- :doc:`tutorial1` - Introduction to PyFEM input files
