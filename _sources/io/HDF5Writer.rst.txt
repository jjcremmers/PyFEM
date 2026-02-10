===========
HDF5Writer
===========

The ``HDF5Writer`` I/O module writes model state and results to an HDF5 file.
In single-file mode it appends each output step to groups named ``cycleN``.
It supports writing nodal displacements, additional nodal fields present in
`globdat.dofs`, custom nodal outputs from `globdat.outputNames`, and optional
per-element outputs.

Overview
--------

Module type: ``HDF5Writer``

- Output file: ``<prefix>.h5`` (prefix from the run `globdat`)
- Mode: single-file with cycle groups (default)
- Methods:
  - ``all``: displacements, extra nodal fields, custom nodal outputs, element outputs
  - ``modes``: modal shapes and eigenvalues (when `globdat.eigenvecs` exists)

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
     - Must be set to ``"HDF5Writer"``

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

Written Data Structure
----------------------

Each cycle group contains:

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Group
     - Datasets (shape → description)
   * - ``elements``
     - ``offsets`` → element offsets; ``connectivity`` → node indices; ``elementIDs``; ``familyIDs``
   * - ``elementGroups``
     - One dataset per group: element IDs
   * - ``nodes``
     - ``coordinates`` → nodal coordinates; ``nodeIDs``
   * - ``nodeGroups``
     - One dataset per group: node IDs
   * - ``nodeData`` (method ``all``)
     - ``displacements`` (N×ndim); extra fields present in `dofTypes` (e.g., ``rx``, ``ry``, ``rz``, ``temp``, ``pres``, ``phase``); custom outputs from `globdat.outputNames`
   * - ``elemData`` (optional)
     - One dataset per element output name when `globdat.elementData` is available
   * - (method ``modes``)
     - ``nodeData/modes`` (nModes×N×ndim); ``eigenvals`` (nModes)

Examples
--------

Minimal output configuration in a ``.pro`` file:

.. code-block:: text

   outputModules = [ "HDF5Writer" ];

   HDF5Writer =
   {
     type     = "HDF5Writer";
     interval = 1;           # write every cycle
   };

See Also
--------
- :doc:`overview`
- Related writers: :doc:`MeshWriter`, :doc:`OutputWriter`, :doc:`GraphWriter`, :doc:`ContourWriter`
