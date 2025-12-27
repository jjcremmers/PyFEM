===========
GraphWriter
===========

The ``GraphWriter`` I/O module writes column-based data per cycle and can
optionally plot the first two columns live to a PNG. It supports data from
``globdat.outputNames``, attributes on ``globdat`` (including arrays indexed
by node/dof), and solver status fields.

Overview
--------

Module type: ``GraphWriter``

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
     - List of column names; for each column, you may define a block with sub-parameters

Per-Column Sub-Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Data source name; defaults to the column name
   * - ``factor``
     - Scale factor applied to the data (default: ``1.0``)
   * - ``node``
     - Node ID, list of node IDs, or node group name when indexing arrays
   * - ``dof``
     - Dof type used to pick the entry from a nodal array (e.g., ``u``, ``v``)

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
     - ``true`` to plot the first two columns live and save ``<prefix>.png`` (default: ``false``)

Examples
--------

- Configure two columns for a force–displacement curve:

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
- :doc:`../io_modules`
- Related writers: ``ContourWriter``, ``OutputWriter``, ``MeshWriter``
