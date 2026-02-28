# I/O Modules

PyFEM provides a comprehensive set of I/O (Input/Output) modules for reading input files, writing simulation results, generating visualizations, and managing data during finite element analyses. These modules enable flexible output of mesh data, nodal results, time histories, and custom visualizations.

## Overview
- **Input modules:** Read problem definitions and restore saved states
- **Output writers:** Export results in various formats (VTK, HDF5, text)
- **Data loggers:** Track specific quantities during analysis (graphs, contours)
- **State management:** Save and restore full analysis states

All output modules share a common configuration pattern and can be combined to produce multiple output formats simultaneously.

## Configuration
Output modules are configured in the `.pro` input file using two components:
1. **Module list:** The `outputModules` property lists which modules to use
2. **Module blocks:** Each module has a configuration block with its parameters

See documentation for configuration examples and details.
