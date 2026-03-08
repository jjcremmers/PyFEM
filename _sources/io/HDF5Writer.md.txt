# HDF5Writer

The `HDF5Writer` I/O module writes model state and results to an HDF5 file. In single-file mode it appends each output step to groups named `cycleN`.

## Overview
- **Module type:** `HDF5Writer`
- Output file: `<prefix>.h5` (prefix from the run `globdat`)
- Mode: single-file with cycle groups (default)
- Methods:
  - `all`: displacements, extra nodal fields, custom nodal outputs, element outputs
  - `modes`: modal shapes and eigenvalues (when `globdat.eigenvecs` exists)

## Parameters
### Mandatory
- `type`: Must be set to `"HDF5Writer"`

### Optional
(See documentation for details)
