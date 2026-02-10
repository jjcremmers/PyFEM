===========
GraphWriter
===========

The ``GraphWriter`` I/O module writes column-based data per cycle and can
optionally plot the first two columns live to a PNG. It supports data from
``globdat.outputNames``, attributes on ``globdat`` (including arrays indexed
by node/dof), and solver status fields.

The output is stored in a multi-column file. By default, the file is stored 
with the extension ``.out``. The first two columns of the output can be shown 
on the screen as a curve during the simulation. Please note that the live 
plotting option only works when PyFEM is used in a Linux environment.

Overview
--------

Module type: ``GraphWriter``

Source: ``pyfem/io/GraphWriter.py``

- Output file: ``<prefix>.out`` by default, or a custom ``filename``.
- Live plot: when ``onScreen = true``, saves a PNG ``<prefix>.png`` each update.
- Columns: configured via a list of column names; each column may have sub-parameters.

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
     - Must be set to ``"GraphWriter"``
   * - ``columns``
     - Array of strings indicating the columns that will be stored. For each 
       column, the type of data, and if needed, the node, degree of freedom 
       and scaling factor needs to be specified.

Per-Column Sub-Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Type of data source. This can be either ``state``, ``velo``, ``fint``, 
       ``stress``, ``lam`` (load factor), etc. Defaults to the column name if not specified.
   * - ``node``
     - Node ID, list of node IDs, or node group name when indexing arrays
   * - ``dof``
     - Degree of freedom. This is most likely ``u`` or ``v``, but can also be 
       ``temp`` for temperature, or ``rx``, ``ry``, ``rz`` for rotations
   * - ``factor``
     - Scale factor applied to the data (default: ``1.0``)

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``filename``
     - Output file name (default: ``<prefix>.out``)
   * - ``onScreen``
     - When set to ``true``, the first two columns will be shown on the screen 
       as a live plot. The default value is ``false``. Note: This only works 
       in a Linux environment.

Examples
--------

Example 1: Load-Displacement Curve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GraphWriter is used to plot the load-displacement curve in the example 
``examples/ch03/NewtonRaphson.pro``. In this case, the external load, the 
load factor ``lambda`` is plotted against the vertical displacement of node 21.

First the graph writer is called in the list of output modules:

.. code-block:: text

   outputModules = ["mesh", "graph"];

where ``mesh`` refers to the VTK output writer, and ``graph`` refers to this 
``GraphWriter``. The following code block in the input file describes the 
input parameters:

.. code-block:: text

   graph =
   {
     type = "GraphWriter";
     onScreen = true;

     columns = [ "disp", "load" ];

     disp = 
     {
       type = "state";
       node = 21;
       dof  = "v";
     };

     load = 
     {
       type = "lam";
     };
   };

In this example, the output file consists of two columns, labeled ``disp`` and 
``load``. The first column contains the displacement of node 21 in the y-direction. 
In the code, the displacement is stored as ``state``, whereas the ``v`` 
degree-of-freedom is the displacement in the vertical direction. In a similar 
fashion, the temperature or a rotation of a specific node can be plotted by 
using the degree-of-freedom names, ``temp``, ``rx``, ``ry`` or ``rz``.

The argument ``onScreen = true;`` indicates that the results are plotted on 
screen during the simulation. Please note that this only works in a Linux environment.

Example 2: Force-Displacement Curve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure two columns for a force-displacement curve:

.. code-block:: text

   outputModules = [ "GraphWriter" ];

   GraphWriter =
   {
     type    = "GraphWriter";
     columns = [ "u", "R" ];

     u = { type = "u"; node = 10; dof = "u"; }       # displacement at node 10
     R = { type = "R"; factor = 1.0; }                # solver reaction or custom output

     onScreen = true;                                   # live plot
     filename = "force_disp.out";
   };

See Also
--------
- :doc:`overview`
- Related writers: :doc:`ContourWriter`, :doc:`OutputWriter`, :doc:`MeshWriter`
