==========
MeshWriter
==========

The ``MeshWriter`` I/O module writes an unstructured VTK file (``.vtu``) per
cycle and a PVD collection file tracking all timesteps. It stores nodes,
elements, nodal dof fields, custom nodal outputs, optional element outputs, and
modal shapes when available.

Overview
--------

Module type: ``MeshWriter``

- Output files: ``<prefix>_t<cycle>.vtu`` and collection ``<prefix>.pvd``.
- Node data: active dof fields and ``globdat.outputNames`` as scalar arrays.
- Element data: optional arrays from ``globdat.elementData``.
- Modes: stores eigenmodes when ``globdat.eigenvecs`` is present.

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
     - Must be set to ``"MeshWriter"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``elementGroup``
     - Element group to write (default: ``"All"``)
   * - ``interval``
     - Output interval in cycles (default: ``1``)
   * - ``extraFields``
     - Additional fields to write (not commonly used; default: empty)
   * - ``beam`` / ``interface``
     - Toggle special handling for these element types (not commonly required)
   * - ``format``
     - ``"binary"`` (default) or ``"ascii"``

Examples
--------

- ``examples/ch03/cantilever8.pro``
- ``examples/elements/sls/sls_cantilever01.pro``

See Also
--------
- :doc:`../io_modules`
- Related writers: ``HDF5Writer``, ``OutputWriter``, ``GraphWriter``
