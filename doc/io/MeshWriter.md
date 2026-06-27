# MeshWriter

The `MeshWriter` I/O module writes an unstructured VTK file (`.vtu`) per cycle and a PVD collection file tracking all timesteps. It stores nodes, elements, nodal dof fields, custom nodal outputs, optional element outputs, and modal shapes when available.

## Overview
- **Module type:** `MeshWriter`
- Output files: `<prefix>_t<cycle>.vtu` and collection `<prefix>.pvd`.
- Node data: active dof fields and `globdat.outputNames` as scalar arrays.
- Element data: optional arrays from `globdat.elementData`.
- Modes: stores eigenmodes when `globdat.eigenvecs` is present.

## Parameters
### Mandatory
- `type`: Must be set to `"MeshWriter"`

### Optional
(See documentation for details)
