# Quickstart

Getting started with PyFEM is straightforward. Follow these steps to install
and run your first simulation.

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/jjcremmers/PyFEM.git
   cd PyFEM
   ```

2. Install PyFEM and its dependencies:

   ```bash
   pip install .
   ```

   This installs the `pyfem` command-line tool and all required packages.

## Verifying the Installation

To verify that the Python environment is working, you can first run a small
standalone example script:

```bash
cd examples/ch02
python PatchTest.py
```

You can then run a full PyFEM model from a `.pro` input file.

## Running a PyFEM Example

After installation, run a basic example to verify everything works:

```bash
cd examples/ch02
pyfem PatchTest8.pro
```

This runs a simple patch test. You should see solver output showing convergence
information and generated result files.

PyFEM includes numerous examples organized by chapter:

```bash
# Navigate to examples
cd examples

# Run a nonlinear truss analysis
cd ch04
pyfem ShallowtrussRiks.pro

# Run a cantilever beam example
cd ../ch03
pyfem cantilever8.pro
```

Each example produces output files (VTK format for visualization, graphs, etc.)
that can be viewed in ParaView or other post-processing tools.

## Understanding a PyFEM Input File

PyFEM input files conventionally end with `.pro`. A typical input file points
to the mesh or data file and defines the active elements and solver settings:

```text
input = "ShallowTrussRiks.dat";
TrussElem  = {
   ...
}
SpringElem = {
   ...
}
solver = {
   ...
}
```

See the user manual and example files for the complete syntax and available
options.

For more detailed installation instructions including virtual environments,
platform-specific notes, and advanced CLI usage, please refer to the
[installation guide](../installation/overview.md).