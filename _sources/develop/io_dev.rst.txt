=========================
Developing I/O Modules
=========================

This guide explains how to implement new input/output modules in PyFEM. I/O
modules handle reading input files, writing results to various formats, and
managing data flow during finite element analysis.

Overview
--------

I/O modules in PyFEM are responsible for:

- Reading problem definitions from input files
- Writing simulation results (VTK, HDF5, text)
- Logging time histories and specific quantities
- Saving and restoring analysis state
- Generating visualizations and plots

All output modules inherit from ``BaseModule`` and implement a ``run`` method
that is called at specified intervals during the analysis.

I/O Module Class Structure
---------------------------

Base Class
~~~~~~~~~~

All I/O modules inherit from ``BaseModule`` located in
``pyfem/util/BaseModule.py``. This provides:

- Automatic property parsing from input files
- Integration with the output manager
- Logging capabilities
- Common utilities

Required Methods
~~~~~~~~~~~~~~~~

An output module must implement:

.. code-block:: python

   class MyOutputModule(BaseModule):
   
       def __init__(self, props, globdat):
           """Initialize output module.
           
           Args:
               props: Properties dictionary from input file
               globdat: Global data container
           """
           BaseModule.__init__(self, props)
           # Initialize module-specific data
           
       def run(self, props, globdat):
           """Write output for current solution step.
           
           Called automatically by the framework at each output interval.
           
           Args:
               props: Properties tree
               globdat: Global data with current solution state
           """
           # Implementation here

Optional Methods
~~~~~~~~~~~~~~~~

For advanced functionality:

.. code-block:: python

   def writeHeader(self):
       """Write file header (called once at initialization)."""
       pass
       
   def finalize(self, globdat):
       """Finalize output (called at end of analysis)."""
       pass

Implementation Examples
-----------------------

Example 1: Text-Based Output Writer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example implements a simple writer that outputs nodal data to a text file:

.. code-block:: python

   # SPDX-License-Identifier: MIT
   # Copyright (c) 2011–2026 Your Name

   from pyfem.util.BaseModule import BaseModule
   import numpy as np

   class SimpleOutputWriter(BaseModule):
       """Write nodal displacements to a text file.
       
       Outputs node ID, coordinates, and displacement components
       in a human-readable format.
       
       Parameters:
           filename: Output file name (default: <prefix>_output.txt)
           interval: Output interval in cycles (default: 1)
           precision: Number of decimal places (default: 6)
       """
   
       def __init__(self, props, globdat):
           """Initialize output writer.
           
           Args:
               props: Properties containing configuration
               globdat: Global data for accessing prefix
           """
           # Default parameters
           self.filename = f"{globdat.prefix}_output.txt"
           self.interval = 1
           self.precision = 6
           
           # Read parameters from input file
           BaseModule.__init__(self, props)
           
           # Open output file and write header
           self.file = open(self.filename, 'w')
           self.writeHeader()
   
       def writeHeader(self):
           """Write file header with column descriptions."""
           self.file.write("# PyFEM Output File\n")
           self.file.write("# Columns: Node | x | y | u | v\n")
           self.file.write("#" + "="*60 + "\n")
   
       def run(self, props, globdat):
           """Write current nodal data.
           
           Outputs data only if current cycle matches the interval.
           
           Args:
               props: Properties tree (not used)
               globdat: Global data with current state
           """
           # Check if we should output this cycle
           if not globdat.solverStatus.cycle % self.interval == 0:
               return
           
           # Write cycle header
           cycle = globdat.solverStatus.cycle
           lam = globdat.lam if hasattr(globdat, 'lam') else 0.0
           
           self.file.write(f"\n# Cycle: {cycle}, Load factor: {lam:.6e}\n")
           
           # Get nodes and DOFs
           nodes = globdat.nodes
           dofs = globdat.dofs
           state = globdat.state
           
           # Loop over all nodes
           for nodeID in nodes.getNodeIDs():
               # Get coordinates
               coords = nodes.getNodeCoords(nodeID)
               
               # Get displacements
               dofIDs = dofs.getForNode(nodeID)
               disps = state[dofIDs]
               
               # Format output
               line = f"{nodeID:6d}"
               
               # Coordinates
               for coord in coords:
                   line += f" {coord:12.{self.precision}e}"
               
               # Displacements
               for disp in disps:
                   line += f" {disp:12.{self.precision}e}"
               
               self.file.write(line + "\n")
           
           # Flush to ensure data is written
           self.file.flush()
   
       def finalize(self, globdat):
           """Close output file.
           
           Called automatically at end of analysis.
           """
           self.file.write("\n# End of output\n")
           self.file.close()

Example 2: Graph Writer for Time Histories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example implements a writer that tracks specific quantities over time:

.. code-block:: python

   from pyfem.util.BaseModule import BaseModule
   import numpy as np
   import matplotlib
   matplotlib.use('Agg')  # Non-interactive backend
   import matplotlib.pyplot as plt

   class TimeHistoryWriter(BaseModule):
       """Write and plot time history data.
       
       Tracks user-specified quantities (displacements, forces, etc.)
       during the analysis and optionally generates plots.
       
       Parameters:
           filename: Output data file (default: <prefix>_history.dat)
           onScreen: Generate PNG plots (default: False)
           columns: List of column names to track
           
           For each column, specify:
               type: Data type ('state', 'fint', 'lam', etc.)
               node: Node ID (for nodal quantities)
               dof: DOF type ('u', 'v', 'w', etc.)
       """
   
       def __init__(self, props, globdat):
           """Initialize time history writer."""
           self.filename = f"{globdat.prefix}_history.dat"
           self.onScreen = False
           self.columns = []
           
           BaseModule.__init__(self, props)
           
           # Storage for time history data
           self.data = {col: [] for col in self.columns}
           self.cycles = []
           
           # Open output file
           self.file = open(self.filename, 'w')
           self.writeFileHeader()
   
       def writeFileHeader(self):
           """Write column headers to data file."""
           header = "# Cycle"
           for col in self.columns:
               header += f"  {col:>15s}"
           self.file.write(header + "\n")
   
       def run(self, props, globdat):
           """Record data for current cycle."""
           cycle = globdat.solverStatus.cycle
           self.cycles.append(cycle)
           
           line = f"{cycle:6d}"
           
           # Extract each column's data
           for colName in self.columns:
               # Get column configuration from props
               colProps = getattr(self, colName)
               
               # Extract value based on type
               value = self.extractValue(colProps, globdat)
               
               # Store and write
               self.data[colName].append(value)
               line += f"  {value:15.6e}"
           
           self.file.write(line + "\n")
           self.file.flush()
           
           # Update plot if requested
           if self.onScreen:
               self.updatePlot()
   
       def extractValue(self, colProps, globdat):
           """Extract a value based on column configuration.
           
           Args:
               colProps: Column properties (type, node, dof)
               globdat: Global data
               
           Returns:
               float: Extracted value
           """
           dataType = colProps.get('type', 'state')
           
           if dataType == 'state':
               # Nodal displacement
               nodeID = colProps['node']
               dofType = colProps['dof']
               dofID = globdat.dofs.getForType(nodeID, dofType)
               return globdat.state[dofID]
               
           elif dataType == 'fint':
               # Internal force
               nodeID = colProps['node']
               dofType = colProps['dof']
               dofID = globdat.dofs.getForType(nodeID, dofType)
               return globdat.fint[dofID]
               
           elif dataType == 'lam':
               # Load factor
               return globdat.lam
               
           elif dataType in globdat.outputNames:
               # Custom output from elements
               data = globdat.getData(dataType, [colProps['node']])
               return data[colProps['node']]
               
           else:
               raise ValueError(f"Unknown data type: {dataType}")
   
       def updatePlot(self):
           """Generate live plot of first two columns."""
           if len(self.columns) < 2:
               return
           
           plt.figure(figsize=(8, 6))
           
           xdata = self.data[self.columns[0]]
           ydata = self.data[self.columns[1]]
           
           plt.plot(xdata, ydata, 'b-o', linewidth=2, markersize=4)
           plt.xlabel(self.columns[0], fontsize=12)
           plt.ylabel(self.columns[1], fontsize=12)
           plt.grid(True, alpha=0.3)
           plt.tight_layout()
           
           # Save plot
           plotname = self.filename.replace('.dat', '.png')
           plt.savefig(plotname, dpi=150)
           plt.close()
   
       def finalize(self, globdat):
           """Close file and generate final plot."""
           self.file.close()
           
           if self.onScreen:
               self.updatePlot()

Example 3: Binary Data Writer (HDF5)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For large-scale simulations, implement binary output using HDF5:

.. code-block:: python

   from pyfem.util.BaseModule import BaseModule
   import h5py
   import numpy as np

   class HDF5Writer(BaseModule):
       """Write simulation data to HDF5 format.
       
       HDF5 provides efficient binary storage for large datasets
       and supports hierarchical data organization.
       
       Parameters:
           interval: Output interval (default: 1)
       """
   
       def __init__(self, props, globdat):
           """Initialize HDF5 writer."""
           self.interval = 1
           BaseModule.__init__(self, props)
           
           # Create HDF5 file
           self.filename = f"{globdat.prefix}.h5"
           self.h5file = h5py.File(self.filename, 'w')
           
           # Write mesh data (written once)
           self.writeMeshData(globdat)
           
           # Counter for time steps
           self.step_count = 0
   
       def writeMeshData(self, globdat):
           """Write mesh topology (written once at initialization)."""
           mesh_group = self.h5file.create_group('mesh')
           
           # Node coordinates
           nodeIDs = globdat.nodes.getNodeIDs()
           nNodes = len(nodeIDs)
           coords = np.zeros((nNodes, globdat.nodes.rank))
           
           for i, nodeID in enumerate(nodeIDs):
               coords[i, :] = globdat.nodes.getNodeCoords(nodeID)
           
           mesh_group.create_dataset('coordinates', data=coords)
           mesh_group.create_dataset('nodeIDs', data=np.array(nodeIDs))
           
           # Element connectivity
           elements = []
           for elemGroup in globdat.elements:
               for elem in elemGroup:
                   elements.append(elem.getNodes())
           
           # Store as variable-length data (different element types)
           mesh_group.create_dataset('connectivity', 
                                    data=np.array(elements, dtype=object))
   
       def run(self, props, globdat):
           """Write solution data for current step."""
           if not globdat.solverStatus.cycle % self.interval == 0:
               return
           
           # Create group for this time step
           step_name = f"step_{self.step_count:06d}"
           step_group = self.h5file.create_group(step_name)
           
           # Write metadata
           step_group.attrs['cycle'] = globdat.solverStatus.cycle
           if hasattr(globdat, 'lam'):
               step_group.attrs['load_factor'] = globdat.lam
           if hasattr(globdat.solverStatus, 'time'):
               step_group.attrs['time'] = globdat.solverStatus.time
           
           # Write displacement field
           step_group.create_dataset('displacements', 
                                    data=globdat.state.copy())
           
           # Write velocity (if available)
           if hasattr(globdat, 'velo'):
               step_group.create_dataset('velocities',
                                        data=globdat.velo.copy())
           
           # Write custom output fields
           for name in globdat.outputNames:
               data = globdat.getData(name, globdat.nodes.getNodeIDs())
               step_group.create_dataset(name, data=np.array(data))
           
           self.step_count += 1
           
           # Flush to disk periodically
           if self.step_count % 10 == 0:
               self.h5file.flush()
   
       def finalize(self, globdat):
           """Close HDF5 file."""
           self.h5file.close()

Example 4: VTK Writer for Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For visualization in ParaView or similar tools:

.. code-block:: python

   from pyfem.util.BaseModule import BaseModule
   import numpy as np

   class VTKWriter(BaseModule):
       """Write VTK files for visualization in ParaView.
       
       Outputs unstructured grid in VTK XML format (.vtu)
       with nodal and element data.
       
       Parameters:
           interval: Output interval (default: 1)
           format: 'binary' or 'ascii' (default: 'binary')
       """
   
       def __init__(self, props, globdat):
           """Initialize VTK writer."""
           self.interval = 1
           self.format = 'binary'
           BaseModule.__init__(self, props)
           
           self.prefix = globdat.prefix
           self.cycle_count = 0
           
           # Create PVD collection file
           self.pvd_file = open(f"{self.prefix}.pvd", 'w')
           self.writePVDHeader()
   
       def writePVDHeader(self):
           """Write ParaView Data (PVD) file header."""
           self.pvd_file.write('<?xml version="1.0"?>\n')
           self.pvd_file.write('<VTKFile type="Collection" version="0.1">\n')
           self.pvd_file.write('  <Collection>\n')
   
       def run(self, props, globdat):
           """Write VTU file for current step."""
           if not globdat.solverStatus.cycle % self.interval == 0:
               return
           
           # Generate filename
           vtu_name = f"{self.prefix}_t{self.cycle_count:06d}.vtu"
           
           # Write VTU file
           self.writeVTU(vtu_name, globdat)
           
           # Add entry to PVD collection
           time = getattr(globdat.solverStatus, 'time', self.cycle_count)
           self.pvd_file.write(
               f'    <DataSet timestep="{time}" file="{vtu_name}"/>\n'
           )
           self.pvd_file.flush()
           
           self.cycle_count += 1
   
       def writeVTU(self, filename, globdat):
           """Write VTK unstructured grid file.
           
           This is a simplified version. Full implementation would
           include proper XML formatting and binary encoding.
           """
           with open(filename, 'w') as f:
               # VTK XML header
               f.write('<?xml version="1.0"?>\n')
               f.write('<VTKFile type="UnstructuredGrid" version="0.1">\n')
               f.write('  <UnstructuredGrid>\n')
               
               # Get mesh data
               nodeIDs = globdat.nodes.getNodeIDs()
               nNodes = len(nodeIDs)
               
               f.write(f'    <Piece NumberOfPoints="{nNodes}" ')
               f.write(f'NumberOfCells="{len(globdat.elements[0])}">\n')
               
               # Point data (nodal fields)
               f.write('      <PointData>\n')
               self.writeNodalData(f, globdat)
               f.write('      </PointData>\n')
               
               # Points (node coordinates)
               f.write('      <Points>\n')
               self.writeNodeCoordinates(f, globdat)
               f.write('      </Points>\n')
               
               # Cells (element connectivity)
               f.write('      <Cells>\n')
               self.writeElementConnectivity(f, globdat)
               f.write('      </Cells>\n')
               
               f.write('    </Piece>\n')
               f.write('  </UnstructuredGrid>\n')
               f.write('</VTKFile>\n')
   
       def writeNodalData(self, f, globdat):
           """Write nodal field data (displacements, etc.)."""
           # Displacement field
           f.write('        <DataArray type="Float64" ')
           f.write('Name="Displacement" NumberOfComponents="3" ')
           f.write('format="ascii">\n')
           
           for nodeID in globdat.nodes.getNodeIDs():
               dofIDs = globdat.dofs.getForNode(nodeID)
               disps = globdat.state[dofIDs]
               
               # Pad to 3D for VTK
               disp_3d = np.zeros(3)
               disp_3d[:len(disps)] = disps
               
               f.write(f"          {disp_3d[0]:.6e} {disp_3d[1]:.6e} ")
               f.write(f"{disp_3d[2]:.6e}\n")
           
           f.write('        </DataArray>\n')
   
       def finalize(self, globdat):
           """Close PVD collection file."""
           self.pvd_file.write('  </Collection>\n')
           self.pvd_file.write('</VTKFile>\n')
           self.pvd_file.close()

Registration and Usage
----------------------

File Location
~~~~~~~~~~~~~

Place your I/O module in:

- ``pyfem/io/MyIOModule.py``

Import in __init__.py
~~~~~~~~~~~~~~~~~~~~~

Add to ``pyfem/io/__init__.py``:

.. code-block:: python

   from .MyIOModule import MyIOModule
   
   __all__ = [
       'MyIOModule',
       # ... other modules
   ]

Using in Input Files
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   outputModules = ["myOutput"];

   myOutput = 
   {
     type = "MyIOModule";
     
     # Module-specific parameters
     filename = "results.txt";
     interval = 5;
   };

Integration with Output Manager
--------------------------------

The output manager automatically:

- Creates instances of all output modules
- Calls ``run()`` at appropriate times
- Handles exceptions
- Calls ``finalize()`` at end of analysis

You don't need to manually register modules with the output manager.

Best Practices
--------------

File Management
~~~~~~~~~~~~~~~

1. **Open files in __init__**: Initialize file handles early
2. **Flush regularly**: Use ``file.flush()`` to ensure data is written
3. **Close in finalize**: Always close files properly
4. **Handle paths**: Use absolute paths or prefix from globdat

Data Formats
~~~~~~~~~~~~

1. **Choose appropriate format**:
   - Text: Human-readable, small datasets
   - Binary (HDF5): Large datasets, efficient
   - VTK: Visualization
2. **Include metadata**: Cycle number, time, load factor
3. **Document format**: Provide header with column descriptions

Error Handling
~~~~~~~~~~~~~~

.. code-block:: python

   def run(self, props, globdat):
       try:
           # Write data
           pass
       except IOError as e:
           logger.error(f"Failed to write output: {e}")
           raise

Performance
~~~~~~~~~~~

1. **Buffer writes**: Don't write every iteration if not needed
2. **Use binary formats**: For large data
3. **Compress data**: HDF5 supports compression
4. **Profile I/O**: Identify bottlenecks

Testing
-------

Unit Tests
~~~~~~~~~~

.. code-block:: python

   import unittest
   import os
   
   class TestMyIOModule(unittest.TestCase):
   
       def test_file_creation(self):
           """Test that output file is created."""
           # Run analysis with module
           # Check file exists
           # Verify content
           pass

Integration Tests
~~~~~~~~~~~~~~~~~

Test with complete analyses to ensure:

- Files are created correctly
- Data is accurate
- Format is valid
- Visualization works (for VTK)

Common Pitfalls
---------------

1. **Not checking interval**: Always check cycle % interval
2. **File handle leaks**: Close files in finalize()
3. **Hardcoded paths**: Use globdat.prefix
4. **Incorrect data types**: Match HDF5/VTK type requirements
5. **Missing error handling**: Catch I/O exceptions

References
----------

For background on finite element output and post-processing:

*"Non-Linear Finite Element Analysis of Solids and Structures"*
by R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
John Wiley & Sons, 2012, ISBN 978-0470666449

External Resources:

- VTK File Formats: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
- HDF5 Documentation: https://docs.h5py.org/
- ParaView Guide: https://www.paraview.org/paraview-guide/

See Also
--------

- :doc:`elements_dev` - Implementing element formulations
- :doc:`materials_dev` - Implementing material models
- :doc:`solvers_dev` - Implementing solution algorithms
- :doc:`../io/overview` - Available I/O modules
