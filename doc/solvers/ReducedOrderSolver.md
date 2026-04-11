# ReducedOrderSolver

The `ReducedOrderSolver` solves the nonlinear equilibrium equations in a
reduced basis generated offline from snapshot data.

## Overview
- **Solver type:** `ReducedOrderSolver`
- **Input data:** HDF5 file containing a `modes` dataset
- **Use case:** online reduced-order solve after snapshot collection and basis construction

## Parameters
- `type`: Must be set to `"ReducedOrderSolver"`
- `modes`: HDF5 file containing the reduced basis
- `modeCount`: Number of basis vectors to keep from the stored modes
- `tol`: Newton tolerance in reduced coordinates
- `iterMax`: Maximum Newton iterations per load step
- `dtime`: Time increment used by the load update
- `loadFunc`: Load scaling function when a load table is not provided

## Example

```text
solver =
{
  type = "ReducedOrderSolver";
  modes = "clamped.h5";
  modeCount = 3;
};
```

## Choosing the Number of Modes

Use the singular values and captured energy reported by `ROMBasisBuilder` to
choose `modeCount`.

- If the first one or two modes already capture most of the energy, a very low-dimensional ROM may be sufficient.
- If the cumulative energy increases slowly, keep more modes.
- Increasing `modeCount` usually improves accuracy at the cost of a larger reduced system.

## Notes
- The number of rows in the stored basis must match the number of global DOFs in the full-order model.
- This solver assumes the basis has already been computed and stored offline.