===========
I/O Modules
===========

PyFEM provides a comprehensive set of I/O (Input/Output) modules for reading 
input files, writing simulation results, generating visualizations, and 
managing data during finite element analyses. These modules enable flexible 
output of mesh data, nodal results, time histories, and custom visualizations.

Overview
--------

The I/O system consists of:

- **Input modules**: Read problem definitions and restore saved states
- **Output writers**: Export results in various formats (VTK, HDF5, text)
- **Data loggers**: Track specific quantities during analysis (graphs, contours)
- **State management**: Save and restore full analysis states

All output modules share a common configuration pattern and can be combined 
to produce multiple output formats simultaneously.

Configuration
-------------

Output modules are configured in the ``.pro`` input file using two components:

1. **Module list**: The ``outputModules`` property lists which modules to use
2. **Module blocks**: Each module has a configuration block with its parameters

Basic Syntax
~~~~~~~~~~~~

.. code-block:: text

   outputModules = ["module1", "module2", "module3"];

   module1 = 
   {
     type = "ModuleType";
     # module-specific parameters
   };

   module2 = 
   {
     type = "AnotherModuleType";
     # module-specific parameters
   };

The name in the ``outputModules`` list corresponds to the configuration block 
name. The ``type`` parameter specifies which I/O module to instantiate.

Complete Example
~~~~~~~~~~~~~~~~

Here's a typical configuration that combines multiple output types:

.. code-block:: text

   outputModules = ["vtk", "graph", "contour"];

   vtk =
   {
     type         = "MeshWriter";
     elementGroup = "All";
     interval     = 1;
   };

   graph =
   {
     type     = "GraphWriter";
     onScreen = true;
     columns  = ["disp", "load"];
     
     disp = { type = "state"; node = 21; dof = "v"; };
     load = { type = "lam"; };
   };

   contour =
   {
     type     = "ContourWriter";
     nodes    = [10, 11, 12, 13, 14, 15];
     interval = 5;
   };

     contour =
   {
     type     = "ContourWriter";
     nodes    = [10, 11, 12, 13, 14, 15];
     interval = 5;
   };

Available Modules
-----------------

Input Modules
~~~~~~~~~~~~~

- :doc:`InputReader` - Reads ``.pro`` files and restores analysis states from dumps
- :doc:`DataDump` - Saves complete analysis state for restart or inspection

Output Writers
~~~~~~~~~~~~~~

- :doc:`MeshWriter` - Writes VTK files (``.vtu``) for visualization in ParaView
- :doc:`HDF5Writer` - Writes HDF5 files for large-scale data storage
- :doc:`OutputWriter` - Writes human-readable text summaries of nodal data

Data Loggers
~~~~~~~~~~~~

- :doc:`GraphWriter` - Tracks time histories and generates plots
- :doc:`ContourWriter` - Outputs data along specified node contours

System Management
~~~~~~~~~~~~~~~~~

- :doc:`OutputManager` - Coordinates all output modules during analysis

Common Parameters
-----------------

Most output modules support these common parameters:

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Module type (mandatory for all modules)
   * - ``interval``
     - Output frequency in analysis cycles (default: 1)
   * - ``filename``
     - Custom output filename (module-dependent default)
   * - ``onScreen``
     - Display output during analysis (where applicable)

Best Practices
--------------

1. **Choose appropriate intervals**: Use larger intervals for expensive output operations
2. **Combine formats**: Use VTK for visualization and text files for time histories
3. **Save binary dumps**: Use DataDump for long analyses to enable restarts
4. **Monitor progress**: Enable GraphWriter with ``onScreen = true`` for feedback

Module Reference
----------------

Click on any module below for detailed documentation, parameters, and examples:

.. toctree::
	:maxdepth: 1
	:titlesonly:

	ContourWriter.rst
	DataDump.rst
	GraphWriter.rst
	HDF5Writer.rst
	InputReader.rst
	MeshWriter.rst
	OutputWriter.rst
