=============
OutputWriter
=============

The ``OutputWriter`` I/O module writes a human-readable summary of all nodes
per cycle to a text file. It uses the global print routine to include node IDs,
coordinates, and active dof values.

Overview
--------

Module type: ``OutputWriter``

- Output file: ``<prefix>_glob.out`` by default, or a custom ``filename``.
- On-screen print: optionally prints the same summary to stdout.

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
     - Must be set to ``"OutputWriter"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``filename``
     - Output file name (default: ``<prefix>_glob.out``)
   * - ``onScreen``
     - ``true`` to also print to stdout (default: ``false``)

Examples
--------

- ``examples/elements/sls/sls_cantilever02.pro``

See Also
--------
- :doc:`../io_modules`
- Related writers: ``MeshWriter``, ``GraphWriter``, ``ContourWriter``
