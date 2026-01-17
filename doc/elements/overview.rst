Elements
========

PyFEM provides a comprehensive library of finite element formulations for 
structural and solid mechanics analysis. Elements define the kinematic 
assumptions, interpolation functions, and integration schemes that transform 
continuum mechanics equations into discrete algebraic systems.

Overview
--------

Element types in PyFEM include:

- **Continuum elements**: 2D, 3D, and axisymmetric solids for bulk materials
- **Structural elements**: Beams, trusses, plates, and shells for slender members
- **Interface elements**: Cohesive zone models for fracture and delamination
- **Special elements**: Springs and other connectors

Each element type is implemented as a Python class that computes element 
stiffness matrices, internal forces, and output quantities. Elements are 
configured in the ``.pro`` input file and associated with specific element 
groups from the mesh.

Configuration in Input Files
-----------------------------

Elements are defined by creating named element groups in the ``.pro`` file. 
Each group specifies the element type and its material properties:

.. code-block:: text

   input = "model.dat";

   ContElem = 
   {
     type = "SmallStrainContinuum";
     
     material = 
     {
       type = "PlaneStress";
       E    = 210.0e3;
       nu   = 0.3;
     };
   };

   solver = 
   {
     type = "NonlinearSolver";
     maxCycle = 50;
   };

In this example:

- ``ContElem`` is the element group name (must match a group in ``model.dat``)
- ``type`` specifies the element formulation
- ``material`` defines the constitutive behavior
- The mesh file ``model.dat`` assigns elements to groups

Multiple Element Groups
~~~~~~~~~~~~~~~~~~~~~~~~

Complex models can combine different element types:

.. code-block:: text

   Solid = 
   {
     type = "SmallStrainContinuum";
     material = { type = "PlaneStress"; E = 210.0e3; nu = 0.3; };
   };

   Beam = 
   {
     type = "Beam";
     material = { type = "Isotropic"; E = 210.0e3; nu = 0.3; };
     A = 0.01;  # Cross-sectional area
     I = 8.33e-6;  # Moment of inertia
   };

   Interface = 
   {
     type = "Interface";
     material = { type = "XuNeedleman"; Tult = 10.0; Gc = 1.0; };
   };

Element Categories
------------------

Continuum Elements
~~~~~~~~~~~~~~~~~~

Solid continuum elements for bulk materials under small or finite strains.
These elements use displacement-based formulations with material behavior 
defined by constitutive laws. Suitable for 2D plane stress/strain, 3D solids, 
and axisymmetric problems.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   finitestrainaxisym.rst
   finitestraincontinuum.rst
   smallstraincontinuum.rst
   smallstrainaxisym.rst   

Shell and Plate Elements
~~~~~~~~~~~~~~~~~~~~~~~~

Elements for thin-walled structures based on Kirchhoff or Reissner-Mindlin 
plate theory. These elements include shear-locking suppression (SLS) variants 
for improved performance with thin structures.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   plate.rst
   sls.rst

Beam and Rod Elements
~~~~~~~~~~~~~~~~~~~~~

One-dimensional structural elements for frames, trusses, and skeletal 
structures. Includes linear and geometrically nonlinear beam formulations 
(Timoshenko, Kirchhoff), truss elements, and spring connectors.

.. toctree::
   :maxdepth: 1
   :titlesonly:

   beamnl.rst
   kirchhoffbeam.rst
   spring.rst
   timoshenkobeam.rst
   truss.rst

Other Elements
~~~~~~~~~~~~~~

Special-purpose elements including interface elements for cohesive zone 
modeling, fracture mechanics, and delamination analysis. These elements use 
traction-separation laws to model progressive failure.

.. toctree::
   :maxdepth: 1
   :titlesonly:
   
   interface.rst
