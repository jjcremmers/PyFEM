# ROMSnapshotWriter

The `ROMSnapshotWriter` module stores full-order state snapshots and the basic
mesh metadata needed for offline reduced-order-model preparation.

## Overview
- **Module type:** `ROMSnapshotWriter`
- **Output format:** HDF5
- **Typical use:** collect snapshot data for POD-based ROM workflows

## Parameters
- `type`: Must be set to `"ROMSnapshotWriter"`
- `filename`: Target HDF5 file used to store the snapshot database

## Behavior
- On the first call it writes:
  - the `state` dataset with the first snapshot
  - element offsets and connectivity
  - nodal coordinates
  - a node-wise displacement DOF map
- On later calls it appends new rows to the `state` dataset

## Example

```text
outputModules = ["vtk", "rom"];

rom =
{
  type = "ROMSnapshotWriter";
  filename = "clamped.h5";
};
```

## Notes
- This module writes snapshots; it does not compute modes itself.
- The resulting HDF5 file is later consumed by `ROMBasisBuilder` and `ReducedOrderSolver`.