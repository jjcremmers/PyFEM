PyFEM Installation Guide
========================

PyFEM can be installed directly from the
`GitHub source <https://github.com/jjcremmers/PyFEM>`_.
Both the **Python API** and the **command-line interface (CLI)** are included.

Requirements
------------

**System Requirements:**

- Python **3.9 or newer**
- ``pip`` (Python package manager)
- Git (for cloning the repository)

**Python Dependencies** (installed automatically):

- ``numpy`` - Numerical computing
- ``scipy`` - Scientific computing
- ``matplotlib`` - Plotting and visualization
- ``meshio`` - Mesh I/O operations
- ``h5py`` - HDF5 file support
- ``PySide6`` - GUI framework
- ``vtk`` - Visualization toolkit

**Recommended: Virtual Environment**

It's recommended to install PyFEM in a virtual environment to avoid conflicts
with other Python packages:

.. code-block:: bash

   # Create virtual environment
   python3 -m venv pyfem-env
   
   # Activate on Linux / macOS
   source pyfem-env/bin/activate
   
   # Activate on Windows PowerShell
   pyfem-env\Scripts\activate
   
   # Activate on Windows Command Prompt
   pyfem-env\Scripts\activate.bat

Installation Steps
------------------

Method 1: Standard Installation (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Clone the repository:**

   .. code-block:: bash

      git clone https://github.com/jjcremmers/PyFEM.git
      cd PyFEM

2. **Install with pip:**

   .. code-block:: bash

      pip install .

   This installs PyFEM and all dependencies, and creates the ``pyfem`` and
   ``pyfem-gui`` command-line executables.

Method 2: Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For developers who want to make local changes and test them immediately:

.. code-block:: bash

   git clone https://github.com/jjcremmers/PyFEM.git
   cd PyFEM
   pip install -e .

The ``-e`` flag installs in "editable" mode, so changes to the source code
are immediately reflected without reinstalling.

Method 3: Direct from GitHub (Advanced)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install directly without cloning:

.. code-block:: bash

   pip install git+https://github.com/jjcremmers/PyFEM.git

Verifying the Installation
---------------------------

After installation, verify that both the CLI and API work correctly.

Checking the CLI
~~~~~~~~~~~~~~~~

The ``pyfem`` command should be available in your terminal:

.. code-block:: bash

   # Check version
   pyfem --help
   
   # Run an example
   cd examples/ch02
   pyfem PatchTest.pro

Expected output includes solver iterations, convergence information, and
generated output files.

Checking the API
~~~~~~~~~~~~~~~~

Test the Python API in an interactive session or script:

.. code-block:: python

   from pyfem import run
   
   # Run a complete analysis
   results = run("examples/ch02/PatchTest.pro")
   
   # Access results
   print(results['globdat'].state)  # Displacement vector

If both tests complete without errors, the installation is successful.

Using PyFEM
-----------

Command-Line Interface (CLI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CLI is the primary way to run PyFEM analyses:

**Basic Usage:**

.. code-block:: bash

   pyfem input_file.pro

**Command-Line Options:**

.. code-block:: bash

   pyfem --help                    # Show help
   pyfem -i input.pro              # Specify input file
   pyfem -d state.dump             # Restart from dump file
   pyfem -p param=value            # Override parameter

**Examples:**

.. code-block:: bash

   # Run a nonlinear analysis
   pyfem examples/ch03/cantilever8.pro
   
   # Restart from saved state
   pyfem -d results_cycle100.dump
   
   # Override material property
   pyfem -i model.pro -p E=210000

Python API
~~~~~~~~~~

The API allows programmatic control of PyFEM analyses from Python scripts.

**Simple Usage - Run to Completion:**

.. code-block:: python

   from pyfem import run
   
   # Run complete analysis
   results = run("input.pro")
   
   # Access global data
   globdat = results['globdat']
   displacements = globdat.state
   
   # Access properties
   props = results['props']

**Advanced Usage - Step-by-Step Control:**

.. code-block:: python

   from pyfem.core.api import PyFEMAPI
   
   # Initialize analysis
   api = PyFEMAPI("input.pro")
   
   # Run step by step
   while api.isActive:
       # Perform one load step
       api.step()
       
       # Access current state
       current_disp = api.globdat.state
       load_factor = api.globdat.lam
       
       # Custom processing or checks
       if load_factor > 5.0:
           print(f"Load factor reached {load_factor}")
   
   # Get final results
   results = api.getResults()
   
   # Clean up
   api.close()

**Loading from Pre-parsed Input:**

.. code-block:: python

   from pyfem.io.InputReader import InputRead
   from pyfem.core.api import PyFEMAPI
   
   # Parse input file
   props, globdat = InputRead("input.pro")
   
   # Modify properties programmatically
   props.solver.tol = 1e-6
   
   # Run with modified properties
   api = PyFEMAPI((props, globdat))
   api.runAll()

**Accessing Results:**

.. code-block:: python

   # After running analysis
   globdat = api.globdat
   
   # Nodal displacements
   displacements = globdat.state
   
   # Nodal coordinates
   node_coords = globdat.nodes.getNodeCoords(nodeID)
   
   # Custom output data
   for name in globdat.outputNames:
       data = globdat.getData(name, node_list)
   
   # Load factor (for nonlinear analysis)
   load_factor = globdat.lam
   
   # Solver status
   cycle = globdat.solverStatus.cycle
   converged = globdat.solverStatus.converged

Updating PyFEM
--------------

To update to the latest version from GitHub:

.. code-block:: bash

   cd PyFEM
   git pull origin main
   pip install --upgrade .

Or if you installed directly from GitHub:

.. code-block:: bash

   pip install --upgrade git+https://github.com/jjcremmers/PyFEM.git

Uninstalling
------------

To remove PyFEM:

.. code-block:: bash

   pip uninstall pyfem

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**1. "pyfem: command not found"**

The executable path may not be in your system PATH. Try:

.. code-block:: bash

   # Find where pip installed the executable
   which pyfem  # Linux/macOS
   where pyfem  # Windows
   
   # Add to PATH or use full path
   ~/.local/bin/pyfem input.pro

**2. Import errors**

If you get import errors when using the API:

.. code-block:: bash

   # Reinstall with dependencies
   pip install --force-reinstall pyfem

**3. VTK or GUI issues**

If visualization or GUI doesn't work:

.. code-block:: bash

   # Install system dependencies (Linux)
   sudo apt-get install libgl1-mesa-glx libxkbcommon-x11-0
   
   # On macOS, ensure XQuartz is installed for GUI

**4. Permission errors during installation**

Use the ``--user`` flag to install in your user directory:

.. code-block:: bash

   pip install --user .

Platform-Specific Notes
-----------------------

Linux
~~~~~

PyFEM works out of the box on most Linux distributions. If using system Python,
you may need to install ``python3-venv``:

.. code-block:: bash

   # Debian/Ubuntu
   sudo apt-get install python3-venv python3-pip
   
   # Fedora/RHEL
   sudo dnf install python3-virtualenv python3-pip

macOS
~~~~~

Install Python 3.9+ using Homebrew if needed:

.. code-block:: bash

   brew install python@3.11
   
For GUI support, install XQuartz:

.. code-block:: bash

   brew install --cask xquartz

Windows
~~~~~~~

1. Install Python 3.9+ from `python.org <https://www.python.org/downloads/>`_
2. Ensure "Add Python to PATH" is checked during installation
3. Use PowerShell or Command Prompt for installation
4. Git for Windows is needed for cloning: `git-scm.com <https://git-scm.com/>`_

Running Examples
----------------

PyFEM includes numerous examples organized by chapter:

.. code-block:: bash

   cd examples
   
   # Chapter 2 - Linear elasticity
   cd ch02
   pyfem PatchTest.pro
   
   # Chapter 3 - Nonlinear analysis
   cd ch03
   pyfem cantilever8.pro
   
   # View results in ParaView
   paraview cantilever8.pvd

Each example directory contains:

- ``.pro`` files: Input files with problem definition
- ``.dat`` files: Mesh files
- Output files generated after running (VTK, text, plots)

Development Setup
-----------------

For contributors and developers:

.. code-block:: bash

   # Clone repository
   git clone https://github.com/jjcremmers/PyFEM.git
   cd PyFEM
   
   # Install in editable mode with dev dependencies
   pip install -e .
   
   # Run tests
   python -m pytest test/
   
   # Check code style
   python -m black pyfem/
   python -m mypy pyfem/

Getting Help
------------

- **Documentation**: https://pyfem.readthedocs.io/
- **GitHub Issues**: https://github.com/jjcremmers/PyFEM/issues
- **Examples**: See the ``examples/`` directory
- **Book**: *"Non-Linear Finite Element Analysis of Solids and Structures"* 
  by de Borst et al., John Wiley & Sons, 2012

Next Steps
----------

After installation:

1. Read the :doc:`../tutorials/quickstart` tutorial
2. Explore examples in the ``examples/`` directory
3. Review the :doc:`../pyfem` for module documentation
4. For development, see the :doc:`../develop/overview`  
