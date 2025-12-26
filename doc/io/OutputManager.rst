==============
OutputManager
==============

``OutputManager`` constructs and runs the configured output modules listed in
``props.outputModules``. It dynamically imports module classes by name and
executes their ``run()`` methods each cycle.

Overview
--------

Manager: ``OutputManager(props, globdat)``

- Configuration: property ``outputModules = [ ... ]`` in the ``.pro`` file.
- For each name in ``outputModules``, an optional block with the same name can
  define parameters, including overriding the module ``type``.

Configuration
-------------

Example ``.pro`` excerpt:

.. code-block:: text

   outputModules = [ "MeshWriter", "OutputWriter", "GraphWriter" ];

   MeshWriter =
   {
     type        = "MeshWriter";
     elementGroup = "All";
     interval    = 1;
   };

   OutputWriter =
   {
     type     = "OutputWriter";
     filename = "summary.out";
   };

   GraphWriter =
   {
     type    = "GraphWriter";
     columns = [ "u", "R" ];
     onScreen = false;
   };

Notes
-----

- If a module block defines ``type``, that class is imported from ``pyfem.io``.
- Errors clearly indicate when a module or class cannot be found.

See Also
--------
- :doc:`../io_modules`
- Individual writer docs: :doc:`MeshWriter <MeshWriter>`, :doc:`OutputWriter <OutputWriter>`, :doc:`GraphWriter <GraphWriter>`, :doc:`ContourWriter <ContourWriter>`, :doc:`HDF5Writer <HDF5Writer>`, :doc:`DataDump <DataDump>`
