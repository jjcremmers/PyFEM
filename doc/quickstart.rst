Getting started with PyFEM is straightforward. Follow these steps to install
and run your first simulation.

Installation
~~~~~~~~~~~~

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/jjcremmers/PyFEM.git
      cd PyFEM

2. Install PyFEM and its dependencies:

   .. code-block:: bash

      pip install .

   This installs the ``pyfem`` command-line tool and all required packages.

Running a Simple Example
~~~~~~~~~~~~~~~~~~~~~~~~

After installation, run a basic example to verify everything works:

.. code-block:: bash

   cd examples/ch02
   pyfem PatchTest.pro

This runs a simple patch test. You should see solver output showing convergence
information and generated result files.

PyFEM includes numerous examples organized by chapter:

.. code-block:: bash

   # Navigate to examples
   cd examples
   
   # Run a nonlinear truss analysis
   cd ch04
   pyfem ShallowTrussRiks.pro
   
   # Run a cantilever beam example
   cd ../ch03
   pyfem cantilever8.pro

Each example produces output files (VTK format for visualization, graphs, etc.)
that can be viewed in ParaView or other post-processing tools.

For more detailed installation instructions including virtual environments,
platform-specific notes, and advanced CLI usage, please refer to the
:doc:`installation guide <installation/overview>`.