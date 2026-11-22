# Quickstart

Getting started with PyFEM is straightforward. Follow these steps to install
and run your first simulation.

## Before you start

Install PyFEM by following the [installation guide](../installation/installation.md).
The examples below assume that the virtual environment is active and that the
`pyfem` command is available.

If you have cloned the repository, run all example commands from its root
directory unless stated otherwise.

## Run your first example

Run a basic patch test from the `examples/ch02` directory:

```bash
cd examples/ch02
pyfem PatchTest8.pro
```

This runs a simple patch test. You should see solver output with convergence
information and generated result files.

To run the standalone Python example in the same directory, use:

```bash
python PatchTest.py
```

The command-line interface also supports:

```bash
pyfem --help                    # Show all options
pyfem --version                 # Show the installed version
pyfem input.pro                 # Run an input file
pyfem -d state.dump             # Restart from a dump file
pyfem input.pro -p E=210000    # Override an input parameter
```

The test shortcuts are documented in the installation guide and can be run
from the repository root:

```bash
pyfem -test
pyfem -coverage
```

PyFEM includes numerous examples organized by chapter:

```bash
# Run a nonlinear truss analysis
cd examples/ch04
pyfem ShallowtrussRiks.pro

# Run a cantilever beam example
cd ../ch03
pyfem cantilever8.pro
```

Each example produces output files, including VTK files for visualization and
other result files. VTK output can be viewed with
[ParaView](https://www.paraview.org/download/):

```bash
paraview path/to/output.pvd
```

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

For platform-specific installation, development setup, testing, coverage, and
troubleshooting, see the [installation guide](../installation/installation.md).
