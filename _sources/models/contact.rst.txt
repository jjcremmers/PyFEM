Contact Mechanics Model
=======================

Overview
--------

The contact model in PyFEM implements a penalty-based contact algorithm for simulating mechanical interactions between deformable bodies and rigid obstacles. The model supports both 2D (disc) and 3D (sphere) contact geometries with a simple penalty method that enforces non-penetration constraints.

The contact model is implemented as a special model type that can be added to any structural analysis to impose contact constraints between the mesh and a moving rigid obstacle.

Theory
------

Penalty Method
~~~~~~~~~~~~~~

The contact algorithm uses the penalty method to enforce contact constraints. When a node penetrates a rigid obstacle, a penalty force is applied:

.. math::

   \mathbf{f}_c = k_p \, g \, \mathbf{n}

where:

- :math:`k_p` is the penalty stiffness parameter
- :math:`g` is the penetration depth (overlap)
- :math:`\mathbf{n}` is the outward normal direction from the obstacle surface

Contact Detection
~~~~~~~~~~~~~~~~~~

For a rigid obstacle (disc in 2D or sphere in 3D) with center :math:`\mathbf{c}` and radius :math:`R`, contact is detected when:

.. math::

   g = R - \|\mathbf{x} - \mathbf{c}\| > 0

where :math:`\mathbf{x}` is the current position of a node. The penetration depth :math:`g` measures how far the node has penetrated into the obstacle.

The outward normal from the obstacle surface is:

.. math::

   \mathbf{n} = \frac{\mathbf{x} - \mathbf{c}}{\|\mathbf{x} - \mathbf{c}\|}

Moving Obstacles
~~~~~~~~~~~~~~~~

The obstacle center can move during the analysis according to a prescribed direction:

.. math::

   \mathbf{c}(t) = \mathbf{c}_0 + \lambda(t) \, \mathbf{d}

where:

- :math:`\mathbf{c}_0` is the initial center position
- :math:`\lambda(t)` is the load factor at time :math:`t`
- :math:`\mathbf{d}` is the prescribed direction vector

Tangent Stiffness
~~~~~~~~~~~~~~~~~

The tangent stiffness contribution from contact is:

.. math::

   \mathbf{K}_c = k_p \, \mathbf{n} \otimes \mathbf{n}

This ensures proper convergence in nonlinear Newton-Raphson iterations.

Input File Configuration
------------------------

Basic Configuration
~~~~~~~~~~~~~~~~~~~

The contact model is defined in the input file using a ``contact`` block:

.. code-block:: python

   contact = {
     type = "disc";           # or "sphere" for 3D
     radius = 1.0;            # obstacle radius
     centre = [5.0, 2.0];     # initial center position
     direction = [0.0, -0.5]; # movement direction
     penalty = 1.0e6;         # penalty stiffness
   };

Parameters
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - ``type``
     - string
     - Contact geometry: ``"disc"`` (2D) or ``"sphere"`` (3D)
   * - ``radius``
     - float
     - Radius of the rigid obstacle
   * - ``centre``
     - array
     - Initial center coordinates ``[x, y]`` or ``[x, y, z]``
   * - ``direction``
     - array
     - Direction of obstacle motion ``[dx, dy]`` or ``[dx, dy, dz]``
   * - ``penalty``
     - float
     - Penalty stiffness parameter :math:`k_p`

Guidelines for Penalty Parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The penalty parameter ``penalty`` should be chosen carefully:

- **Too small**: Excessive penetration, inaccurate contact enforcement
- **Too large**: Numerical ill-conditioning, convergence difficulties
- **Recommended**: :math:`k_p \approx 10^2 - 10^4 \times E` where :math:`E` is the material stiffness

Examples
--------

Example 1: 2D Disc Contact
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 2D finite strain analysis with a moving disc obstacle:

.. code-block:: python

   input = "contact_test.dat";

   ContElem = {
     type = "FiniteStrainContinuum";
     material = {
       type = "PlaneStress";
       E = 1.0e6;
       nu = 0.25;
     };
   };

   contact = {
     type = "disc";
     radius = 1.0;
     centre = [5.0, 2.0];
     direction = [0.0, -0.5];
     penalty = 1.0e6;
   };

   solver = {
     type = "NonlinearSolver";
     maxCycle = 20;
     dtime = 0.1;
   };

   outputModules = ["vtk"];

   vtk = {
     type = "MeshWriter";
   };

This example simulates a deformable structure being pressed by a rigid disc moving downward. The disc center starts at (5.0, 2.0) and moves in direction (0.0, -0.5) scaled by the load factor.

**Run the example:**

.. code-block:: bash

   cd examples/contact
   pyfem contact_test.pro
   paraview contact_test.pvd

Example 2: 3D Sphere Contact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 3D small strain analysis with a spherical obstacle:

.. code-block:: python

   input = "contact_test3D_01.dat";

   ContElem = {
     type = "SmallStrainContinuum";
     material = {
       type = "Isotropic";
       E = 1.0e6;
       nu = 0.25;
     };
   };

   contact = {
     type = "sphere";
     radius = 1.0;
     centre = [0.0, 0.0, 1.0];
     direction = [0.0, 0.0, -0.05];
     penalty = 1.0e6;
   };

   solver = {
     type = "NonlinearSolver";
     maxCycle = 20;
     dtime = 0.1;
   };

   outputModules = ["vtk"];

   vtk = {
     type = "MeshWriter";
   };

This example demonstrates 3D contact with a sphere moving in the negative z-direction.

**Run the example:**

.. code-block:: bash

   cd examples/contact
   pyfem contact_test3D_01.pro
   paraview contact_test3D_01.pvd

Example 3: Contact with Output Monitoring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Contact analysis with displacement tracking:

.. code-block:: python

   input = "contact_test02.dat";

   ContElem = {
     type = "FiniteStrainContinuum";
     material = {
       type = "PlaneStress";
       E = 1.0e6;
       nu = 0.25;
     };
   };

   contact = {
     type = "disc";
     radius = 1.0;
     centre = [8.0, 1.4];
     direction = [0.0, -0.5];
     penalty = 1.0e6;
   };

   solver = {
     type = "NonlinearSolver";
     maxCycle = 40;
     dtime = 0.2;
   };

   outputModules = ["vtk", "GraphWriter"];

   vtk = {
     type = "MeshWriter";
   };

   GraphWriter = {
     onScreen = true;
     columns = ["time", "disp"];
     disp = {
       type = "state";
       node = 100;
       dof  = "v";
     };
   };

This example tracks the vertical displacement of a specific node during contact.

Usage Guidelines
----------------

Mesh Requirements
~~~~~~~~~~~~~~~~~

- **Element compatibility**: Works with any continuum element formulation
- **Refinement**: Use adequate mesh refinement near expected contact zones
- **Node distribution**: Ensure sufficient nodes in contact region for accurate penalty enforcement

Solver Configuration
~~~~~~~~~~~~~~~~~~~~

Contact problems require nonlinear solvers:

.. code-block:: python

   solver = {
     type = "NonlinearSolver";
     maxCycle = 20;          # Sufficient iterations for convergence
     dtime = 0.1;            # Time step size
     tol = 1.0e-3;           # Convergence tolerance (optional)
   };

Best Practices
~~~~~~~~~~~~~~

1. **Start with small time steps**: Contact can cause sudden stiffness changes
2. **Monitor convergence**: Check Newton-Raphson iteration counts
3. **Tune penalty parameter**: Balance accuracy and conditioning
4. **Gradual loading**: Increase contact depth gradually through load factor
5. **Check penetration**: Verify that penetration depths remain small

Limitations
-----------

Current Implementation
~~~~~~~~~~~~~~~~~~~~~~

The current contact model has the following limitations:

- **Rigid obstacles only**: No deformable-deformable contact
- **Simple geometries**: Only discs (2D) and spheres (3D)
- **No friction**: Frictionless contact only
- **Penalty method**: Small penetration occurs (not a true constraint)
- **Single obstacle**: One contact object per analysis

Future Extensions
~~~~~~~~~~~~~~~~~

Potential enhancements include:

- Multiple contact objects
- Arbitrary contact surface geometries
- Friction models (Coulomb friction)
- Contact between deformable bodies
- Augmented Lagrangian or Lagrange multiplier methods

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Problem**: Excessive penetration

- **Solution**: Increase ``penalty`` parameter
- **Check**: Verify obstacle position and movement direction

**Problem**: Convergence difficulties

- **Solution**: Reduce ``penalty`` parameter or decrease time step
- **Check**: Ensure mesh quality near contact zone

**Problem**: No contact detected

- **Solution**: Verify ``centre``, ``radius``, and mesh geometry
- **Check**: Obstacle should intersect the deformed mesh

**Problem**: Solver divergence

- **Solution**: Reduce ``dtime`` or increase ``maxCycle``
- **Check**: Initial configuration should have minimal penetration

References
----------

For theoretical background on contact mechanics and penalty methods, see:

- **Chapter 13** of *Non-Linear Finite Element Analysis of Solids and Structures* (de Borst et al., 2012)
- P. Wriggers, *Computational Contact Mechanics*, 2nd Edition, Springer, 2006
- T.A. Laursen, *Computational Contact and Impact Mechanics*, Springer, 2002

See Also
--------

- :doc:`overview` - Models overview
- :doc:`rve` - RVE model for multi-scale analysis
- :doc:`../elements/overview` - Element formulations
- :doc:`../solvers/NonlinearSolver` - Nonlinear solution algorithms
- :doc:`../io/MeshWriter` - VTK output for visualization

**Example files:**

- ``examples/contact/contact_test.pro`` - 2D disc contact
- ``examples/contact/contact_test02.pro`` - 2D with output monitoring
- ``examples/contact/contact_test3D_01.pro`` - 3D sphere contact
