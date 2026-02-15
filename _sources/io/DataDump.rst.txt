=========
DataDump
=========

The ``DataDump`` I/O module serializes the entire analysis state (``props`` and
``globdat``) to a pickle file for later restart or inspection.

Overview
--------

Module type: ``DataDump``

- Output file: ``<prefix>_<cycle>.dump`` or ``<prefix>.dump`` when ``lastOnly``.
- Contents: both the properties object and global data.

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
     - Must be set to ``"DataDump"``

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
   * - ``lastOnly``
     - ``true`` to always overwrite a single dump file (default: ``false``)

Examples
--------

- Save dumps every cycle:

.. code-block:: text

   outputModules = [ "DataDump" ];

   DataDump =
   {
     type     = "DataDump";
     interval = 1;
     lastOnly = false;
   };

See Also
--------
- :doc:`overview`
- :doc:`InputReader` for reading dumps
