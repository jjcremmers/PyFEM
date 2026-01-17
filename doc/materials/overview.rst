Materials
=========

PyFEM provides a range of material models for linear and nonlinear analysis,
including elastic, elastoplastic, and cohesive zone models. Material models
define the constitutive behavior that relates stress to strain (or traction
to displacement jump for cohesive elements).

Configuration in Input Files
-----------------------------

Material models are configured within element definitions in the ``.pro`` 
input file. Each element group can have its own material definition, or 
materials can be defined globally and referenced by name.

Basic Syntax
~~~~~~~~~~~~

The material block is nested inside an element definition:

.. code-block:: text

   ElementGroup = 
   {
     type = "SmallStrainContinuum";
     
     material = 
     {
       type = "PlaneStress";
       E    = 100.0;
       nu   = 0.3;
     };
   };

The ``type`` parameter specifies which material model to use, and subsequent
parameters define the material properties (e.g., Young's modulus ``E``, 
Poisson's ratio ``nu``).

Material Types and State Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PyFEM uses a hierarchical material system:

1. **State models** (e.g., PlaneStress, PlaneStrain) - Define the stress state
2. **Constitutive models** (e.g., HookesLaw, VonMises) - Define material behavior

For simple elastic materials, the state model can be specified directly:

.. code-block:: text

   material = 
   {
     type = "PlaneStress";  # State model only
     E    = 100.0;
     nu   = 0.3;
   };

For complex behavior, nest a constitutive model inside a state model:

.. code-block:: text

   material = 
   {
     type = "PlaneStrain";
     E    = 210.0e3;
     nu   = 0.3;
     
     model = 
     {
       type  = "IsotropicHardeningPlasticity";
       sY    = 250.0;
       hard  = 1000.0;
     };
   };

Named Materials
~~~~~~~~~~~~~~~

For complex models or when using the same material in multiple elements, 
define materials by name and reference them:

.. code-block:: text

   Steel = 
   {
     type = "PlaneStress";
     E    = 210.0e3;
     nu   = 0.3;
   };

   Aluminum = 
   {
     type = "PlaneStress";
     E    = 70.0e3;
     nu   = 0.33;
   };

   Element1 = 
   {
     type     = "SmallStrainContinuum";
     material = "Steel";
   };

   Element2 = 
   {
     type     = "SmallStrainContinuum";
     material = "Aluminum";
   };

Multi-Material Support
~~~~~~~~~~~~~~~~~~~~~~

Some elements support spatially varying materials using ``MultiMaterial``:

.. code-block:: text

   material = 
   {
     type      = "MultiMaterial";
     materials = ["Steel", "Aluminum"];
   };

The material assignment is typically based on element IDs or integration 
point locations defined in the mesh file.

Material Categories
-------------------

Material Categories
-------------------

Constitutive Relations for Continua
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These materials are used with continuum elements (2D, 3D, and axisymmetric):

- Linear elastic models (Hooke's law, isotropic, transverse isotropic)
- Stress state wrappers (plane stress, plane strain)
- Elastoplastic models (isotropic hardening, kinematic hardening)
- Specialized models (sandwich cores, multi-material)

.. toctree::
	:maxdepth: 1
	:titlesonly:

	hookeslaw.rst
    isotropic.rst
	isotropichardeningplasticity.rst
    isokinematichardening.rst
	multimaterial.rst
	planestrain.rst
	planestress.rst
	sandwichcore.rst
	transverseisotropic.rst



Cohesive Constitutive Relations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These materials model interface behavior and fracture mechanics:

- Cohesive zone models for delamination and crack propagation
- Traction-separation laws for interface elements
- Mode I and mixed-mode fracture models

.. toctree::
	:maxdepth: 1
	:titlesonly:

    dummy.rst
    powerlawmodei.rst
	thoulessmodei.rst
	xuneedleman.rst    

Failure Criteria
~~~~~~~~~~~~~~~~

These models define yield surfaces and failure conditions:

- Von Mises plasticity
- Other stress-based failure criteria

.. toctree::
	:maxdepth: 1
	:titlesonly:

	vonmises.rst    