# ContourWriter

The `ContourWriter` I/O module writes a plain-text table of nodal data for a selected set of nodes at a given interval. Each file is named `<prefix>-contour-<k>.out` where `k` increments per write.

## Overview
- **Module type:** `ContourWriter`
- Writes node ID, coordinates, all active dof values, and custom outputs in `globdat.outputNames`.
- One row per node; header includes columns for coordinates and outputs.

## Parameters
### Mandatory
- `type`: Must be set to `"ContourWriter"`
- `nodes`: List of node IDs or a node group name selecting the nodes to output

### Optional
(See documentation for details)
