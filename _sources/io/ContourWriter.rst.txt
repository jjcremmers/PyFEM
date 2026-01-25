=============
ContourWriter
=============

The ``ContourWriter`` I/O module writes a plain-text table of nodal data for a
selected set of nodes at a given interval. Each file is named
``<prefix>-contour-<k>.out`` where ``k`` increments per write.

Overview
--------

Module type: ``ContourWriter``

- Writes node ID, coordinates, all active dof values, and custom outputs in
  ``globdat.outputNames``.
- One row per node; header includes columns for coordinates and outputs.

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
     - Must be set to ``"ContourWriter"``
   * - ``nodes``
     - List of node IDs or a node group name selecting the nodes to output

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``interval``
     - Output interval in cycles (default: ``1``)

Examples
--------

Basic Usage
~~~~~~~~~~~

Here's a complete example showing how to configure the ``ContourWriter`` to
output data along a contour of nodes during a nonlinear analysis:

.. code-block:: python

   input = "model.dat";

   solver = 
   {
     type = "NonlinearSolver";
     maxCycle = 100;
   };

   outputModules = ["contour"];

   contour =
   {
     type = "ContourWriter";
     
     # Specify nodes to track (can be a list or node group)
     nodes = [10, 11, 12, 13, 14, 15];
     
     # Output every 5 cycles
     interval = 5;
   };

This configuration will create output files named ``<prefix>-contour-0.out``,
``<prefix>-contour-1.out``, etc., at every 5th cycle. Each file contains a
table with the following columns:

- Node ID
- Nodal coordinates (x, y, and optionally z)
- All active degree of freedom values (e.g., u, v, w)
- Any custom output variables defined in ``globdat.outputNames``

Using Node Groups
~~~~~~~~~~~~~~~~~

Instead of listing individual nodes, you can reference a node group defined
in your input file:

.. code-block:: python

   contour =
   {
     type = "ContourWriter";
     
     # Reference a node group from the mesh
     nodes = "crackPath";
     
     interval = 1;
   };

This is particularly useful for tracking complex contours or when the node
IDs are not known in advance.

Output Format
~~~~~~~~~~~~~

The output files are plain-text tables with the following format::

   #Node  x-coor     y-coor     u          v          stress     strain    
      10  0.000e+00  1.000e+00  1.234e-03  5.678e-04  2.345e+01  1.234e-03
      11  1.000e-01  1.000e+00  1.456e-03  5.890e-04  2.456e+01  1.345e-03
      12  2.000e-01  1.000e+00  1.678e-03  6.012e-04  2.567e+01  1.456e-03
      ...

This format can be easily imported into post-processing tools or plotting
libraries like matplotlib, gnuplot, or spreadsheet applications.

Additional Examples
~~~~~~~~~~~~~~~~~~~

- ``examples/elements/interface/PeelTest.pro``

See Also
--------
- :doc:`overview`
- Related writers: :doc:`GraphWriter`, :doc:`OutputWriter`, :doc:`MeshWriter`
