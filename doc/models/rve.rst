===================================
RVE (Representative Volume Element)
===================================

The ``RVE`` model implements periodic boundary conditions for representative
volume element analyses, commonly used in computational homogenization and
multi-scale modeling of heterogeneous materials.

Overview
--------

Module type: ``RVE``

Representative Volume Elements (RVEs) are used to determine effective
(homogenized) material properties of heterogeneous materials by analyzing a
small representative sample with periodic microstructure. The RVE model
enforces:

- Periodic displacement fields on opposite boundaries
- Prescribed macroscopic strain states
- Proper constraint coupling for homogenization

This enables the computation of effective material properties from microscale
simulations.

Theoretical Background
----------------------

In computational homogenization, the macroscopic stress-strain relationship is
derived from the microscopic response of an RVE subjected to periodic boundary
conditions. The key principle is that:

1. **Macroscopic strain** :math:`\bar{\boldsymbol{\varepsilon}}` is prescribed
2. **Microscopic displacement** field :math:`\mathbf{u}(\mathbf{x})` is periodic
3. **Macroscopic stress** :math:`\bar{\boldsymbol{\sigma}}` is volume-averaged

For a rectangular RVE, periodic boundary conditions require:

.. math::

   \mathbf{u}(\mathbf{x} + \mathbf{L}_i) - \mathbf{u}(\mathbf{x}) = 
   \bar{\boldsymbol{\varepsilon}} \cdot \mathbf{L}_i

where :math:`\mathbf{L}_i` are the RVE edge vectors.

Parameters
----------

Mandatory Parameters
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Must be set to ``"Periodic"`` or ``"Prescribed"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

Currently, the RVE model uses hardcoded strain values. Future versions will
support strain input parameters.

Boundary Condition Types
-------------------------

The RVE model supports two types of boundary conditions:

Periodic Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Type:** ``"Periodic"``

This is the standard approach for computational homogenization. The displacement
field is decomposed as:

.. math::

   \mathbf{u}(\mathbf{x}) = \bar{\boldsymbol{\varepsilon}} \cdot \mathbf{x} + 
   \mathbf{u}^*(\mathbf{x})

where :math:`\mathbf{u}^*` is the periodic fluctuation field satisfying:

.. math::

   \mathbf{u}^*(\mathbf{x}_{\text{right}}) = \mathbf{u}^*(\mathbf{x}_{\text{left}})
   
   \mathbf{u}^*(\mathbf{x}_{\text{top}}) = \mathbf{u}^*(\mathbf{x}_{\text{bottom}})

**Implementation:**

1. Bottom-left corner node is fixed: :math:`\mathbf{u}_1 = \mathbf{0}`
2. Other corner nodes prescribed based on strain and RVE dimensions
3. Interior boundary nodes coupled with linear constraints

**Advantages:**

- Allows periodic fluctuations in the displacement field
- Provides most accurate homogenization results
- Suitable for unit cells with periodic microstructure

Prescribed Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Type:** ``"Prescribed"``

This simpler approach prescribes displacements at all boundary nodes according
to a linear displacement field:

.. math::

   u_x = \bar{\varepsilon}_{xx} \, x + \frac{1}{2} \bar{\gamma}_{xy} \, y
   
   u_y = \bar{\varepsilon}_{yy} \, y + \frac{1}{2} \bar{\gamma}_{xy} \, x

**Advantages:**

- Simpler to implement and understand
- No constraint coupling required
- Suitable for validation and testing

**Disadvantages:**

- Does not allow periodic fluctuations
- Less accurate for strongly heterogeneous materials
- Not recommended for computational homogenization

Mesh Requirements
-----------------

For the RVE model to work correctly, the mesh must satisfy:

Node Groups
~~~~~~~~~~~

Four node groups must be defined in the mesh file:

- ``Left``: Nodes on the left boundary
- ``Right``: Nodes on the right boundary  
- ``Bottom``: Nodes on the bottom boundary
- ``Top``: Nodes on the top boundary

**Important:** Corner nodes should be included in both relevant groups 
(e.g., bottom-left corner in both ``Left`` and ``Bottom``).

Geometric Requirements
~~~~~~~~~~~~~~~~~~~~~~

- The RVE must be rectangular
- Opposite boundaries must have matching node positions
- Node spacing should be consistent on opposite sides
- The model automatically validates these requirements at initialization

Examples
--------

Example 1: Periodic RVE Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here's a complete example for a periodic RVE analysis:

.. code-block:: text

   input = "rve_mesh.dat";

   ContElem =
   {
     type = "SmallStrainContinuum";

     material =
     {
       type = "PlaneStress";
       E    = 1.0e6;
       nu   = 0.3;
     };
   };

   RVE =
   {
     type = "Periodic";
   };

   solver =
   {
     type = "LinearSolver";
   };

   outputModules = ["vtk"];

   vtk =
   {
     type = "MeshWriter";
   };

This configuration:

1. Loads a mesh with node groups for boundaries
2. Uses periodic boundary conditions
3. Solves for the microscopic displacement field
4. Outputs results to VTK format for visualization

Example 2: Prescribed RVE Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For simpler prescribed boundary conditions:

.. code-block:: text

   input = "rve_mesh.dat";

   ContElem =
   {
     type = "SmallStrainContinuum";

     material =
     {
       type = "PlaneStrain";
       E    = 210.0e3;
       nu   = 0.3;
     };
   };

   RVE =
   {
     type = "Prescribed";
   };

   solver =
   {
     type = "LinearSolver";
   };

   outputModules = ["vtk", "graph"];

   vtk =
   {
     type = "MeshWriter";
   };

Mesh File Format
----------------

The mesh file must define node groups for boundaries:

.. code-block:: text

   <NodeGroup name="Left">
     1 5 9 13
   </NodeGroup>

   <NodeGroup name="Right">
     4 8 12 16
   </NodeGroup>

   <NodeGroup name="Bottom">
     1 2 3 4
   </NodeGroup>

   <NodeGroup name="Top">
     13 14 15 16
   </NodeGroup>

The model automatically sorts nodes and validates the geometry.

Computing Effective Properties
-------------------------------

After running an RVE analysis, effective material properties can be computed
from the homogenized stress-strain relationship.

Macroscopic Stress
~~~~~~~~~~~~~~~~~~

The macroscopic stress is computed as the volume average:

.. math::

   \\bar{\\boldsymbol{\\sigma}} = \\frac{1}{V} \\int_V \\boldsymbol{\\sigma}(\\mathbf{x}) \\, dV

In practice, this is computed by summing element contributions:

.. math::

   \\bar{\\boldsymbol{\\sigma}} = \\frac{1}{V} \\sum_e \\boldsymbol{\\sigma}_e \\, V_e

Effective Stiffness Tensor
~~~~~~~~~~~~~~~~~~~~~~~~~~~

By applying multiple strain states and recording the resulting stresses, the
effective stiffness tensor can be determined:

.. math::

   \\bar{\\boldsymbol{\\sigma}} = \\mathbf{C}^{\\text{eff}} : \\bar{\\boldsymbol{\\varepsilon}}

For 2D plane stress:

1. Apply :math:`\\bar{\\varepsilon}_{xx} = 1, \\bar{\\varepsilon}_{yy} = 0, \\bar{\\gamma}_{xy} = 0`
2. Record :math:`\\bar{\\sigma}_{xx}, \\bar{\\sigma}_{yy}, \\bar{\\sigma}_{xy}`
3. Repeat for other strain states
4. Assemble effective stiffness matrix

Output and Post-Processing
---------------------------

The RVE analysis generates:

Displacement Field
~~~~~~~~~~~~~~~~~~

The microscopic displacement field shows the deformation pattern within the RVE:

- Visualize in ParaView using VTK output
- Check periodicity by examining opposite boundaries
- Identify high-strain regions

Stress Field
~~~~~~~~~~~~

The stress distribution reveals:

- Stress concentrations due to microstructure
- Load transfer mechanisms
- Critical regions for failure analysis

Volume Averaging
~~~~~~~~~~~~~~~~

Compute volume-averaged quantities using:

.. code-block:: python

   from pyfem import run
   
   # Run RVE analysis
   results = run("rve_model.pro")
   globdat = results['globdat']
   
   # Access element stresses (requires custom processing)
   # Compute volume average
   total_volume = sum(element_volumes)
   avg_stress = sum(stress * volume for stress, volume in zip(stresses, volumes)) / total_volume

Limitations and Considerations
-------------------------------

Current Limitations
~~~~~~~~~~~~~~~~~~~

1. **Hardcoded strain**: Strain components are currently hardcoded in the source
2. **2D only**: Only 2D plane stress/strain is supported
3. **Rectangular RVE**: Only rectangular domains are supported
4. **Linear analysis**: Currently works with linear solvers only

Future Extensions
~~~~~~~~~~~~~~~~~

Planned enhancements include:

- Strain input parameters in ``.pro`` file
- 3D RVE support
- Nonlinear RVE analysis
- Multiple strain cases in one run
- Automated effective property computation

Best Practices
--------------

1. **Mesh quality**: Use regular, well-structured meshes
2. **Size effects**: Verify RVE size is sufficient (convergence study)
3. **Periodicity**: Ensure microstructure is truly periodic at boundaries
4. **Validation**: Compare with analytical solutions for simple cases
5. **Multiple strain states**: Apply all independent strain components

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**"Different number of nodes" warning**

- Check that mesh is truly rectangular
- Verify node groups contain correct nodes
- Ensure corner nodes are in both relevant groups

**"Coordinates do not match" warning**

- Opposite boundaries must have matching node positions
- Check mesh generation and ensure proper pairing
- Use consistent discretization on opposite sides

**Singular stiffness matrix**

- Verify boundary conditions are properly applied
- Check that corner nodes are correctly constrained
- Ensure node groups are properly defined

**Non-physical results**

- Verify material properties are reasonable
- Check strain magnitude is appropriate
- Inspect stress distribution for anomalies

Additional Examples
-------------------

Full example files are provided in:

- ``examples/models/rve/rve_test01.pro``
- ``examples/models/rve/rve_test01.dat``

References
----------

For theoretical background on computational homogenization:

- Kouznetsova, V., Brekelmans, W.A.M., and Baaijens, F.P.T. (2001). 
  "An approach to micro-macro modeling of heterogeneous materials." 
  *Computational Mechanics*, 27(1), 37-48.

- Geers, M.G.D., Kouznetsova, V.G., and Brekelmans, W.A.M. (2010). 
  "Multi-scale computational homogenization: Trends and challenges." 
  *Journal of Computational and Applied Mathematics*, 234(7), 2175-2182.

- Miehe, C., Schröder, J., and Schotte, J. (1999). 
  "Computational homogenization analysis in finite plasticity." 
  *Computer Methods in Applied Mechanics and Engineering*, 171(3-4), 387-418.

See Also
--------

- :doc:`overview` - Overview of all models
- :doc:`../elements/overview` - Element formulations
- :doc:`../materials/overview` - Material models
- :doc:`../solvers/overview` - Solution algorithms
