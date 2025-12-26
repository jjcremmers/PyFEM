======================
SmallStrainContinuum
======================

The ``SmallStrainContinuum`` element is a standard displacement-based finite 
element for 2D and 3D small strain continuum analysis. This element formulation 
is suitable for linear and non-linear material behavior under the assumption of 
small strains and small displacements. It supports various material models 
including elastic, plastic, and damage constitutive laws.

The element automatically adapts to 2D or 3D analysis based on the problem 
dimension (``rank``). In 2D, it provides displacement degrees of freedom ``u`` 
and ``v`` for the x and y directions. In 3D, it includes an additional degree 
of freedom ``w`` for the z direction.

--------
Overview
--------

Element type: ``SmallStrainContinuum``

The element implements:

- **2D formulation**: Plane stress or plane strain analysis with triangular 
  (3-node, 6-node) or quadrilateral (4-node, 8-node, 9-node) elements
- **3D formulation**: Full 3D analysis with tetrahedral or hexahedral elements
- **Material nonlinearity**: Compatible with elastic, plastic, damage, and 
  other constitutive models
- **Numerical integration**: Gauss quadrature based on element type

The element computes the tangent stiffness matrix, internal force vector, 
mass matrix (for dynamic analysis), and dissipation (for inelastic materials).

----------
Parameters
----------

Mandatory parameters
~~~~~~~~~~~~~~~~~~~~

``type``
  Must be set to ``"SmallStrainContinuum"``

``material``
  A material block defining the constitutive behavior. The material block must 
  include a ``type`` parameter specifying the material model.

  Common material types include:

  - ``"PlaneStress"``: For 2D plane stress analysis
  - ``"PlaneStrain"``: For 2D plane strain analysis  
  - ``"Isotropic"``: For 3D isotropic elasticity
  - Other material models (e.g., plasticity, damage) as available

  Within the material block, material-specific parameters must be provided, 
  such as:

  - ``E``: Young's modulus
  - ``nu``: Poisson's ratio
  - Additional parameters depending on the material type (e.g., yield stress, 
    hardening parameters)

Optional parameters
~~~~~~~~~~~~~~~~~~~

The ``SmallStrainContinuum`` element itself does not require additional optional 
parameters beyond the mandatory material definition. However, specific material 
models may have optional parameters (refer to the material documentation).

--------
Examples
--------

Example 1: 2D Plane Stress Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 2D continuum element with plane stress material::

  ContElem =
  {
    type = "SmallStrainContinuum";

    material =
    {
      type = "PlaneStress";
      E    = 1.e6;
      nu   = 0.25;
    };
  };

This configuration is used in the patch test example: 
``examples/ch02/PatchTest4.pro``

Example 2: 3D Isotropic Elasticity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 3D continuum element with isotropic material::

  ContElem =
  {
    type = "SmallStrainContinuum";

    material =
    {
      type = "Isotropic";
      E    = 1.e6;
      nu   = 0.25;
    };
  };

This configuration is used in the 3D patch test example: 
``examples/ch02/PatchTest8_3D.pro``

Example 3: Nonlinear Material Behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The element can be combined with more advanced material models for nonlinear 
analysis. For damage mechanics::

  ContElem =
  {
    type = "SmallStrainContinuum";

    material =
    {
      type = "ContinuumDamage";
      E    = 30000.;
      nu   = 0.2;
      fc   = 20.;
      Gf   = 0.1;
    };
  };

This configuration is used in: ``examples/ch06/ContDamExample.pro``

------------------
Additional Examples
------------------

The ``SmallStrainContinuum`` element is used extensively throughout the 
examples directory:

- **Patch tests**: ``examples/ch02/PatchTest*.pro``
- **Interface elements**: ``examples/ch13/PeelTest*.pro``
- **Material models**: ``examples/ch06/ContDamExample.pro``, 
  ``examples/materials/plasticity/dogbone.pro``
- **Gmsh integration**: ``examples/gmsh/twist.pro``, 
  ``examples/gmsh/two_fibres.pro``
- **Contact mechanics**: ``examples/contact/contact_test3D_01.pro``

See Also
--------

- :doc:`materials` - Available material models
- :doc:`tutorial1` - Introduction to PyFEM input files
- :doc:`elements` - Overview of all element types
