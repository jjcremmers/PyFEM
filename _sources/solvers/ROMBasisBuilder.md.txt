# ROMBasisBuilder

The `ROMBasisBuilder` solver constructs a reduced basis from previously
collected full-order snapshots stored in an HDF5 file.

## Overview
- **Solver type:** `ROMBasisBuilder`
- **Current method:** POD
- **Input data:** snapshot matrix stored in the `state` dataset of an HDF5 file
- **Output data:** `modes` and `eigenvals` datasets written to the same HDF5 file

## Parameters
- `type`: Must be set to `"ROMBasisBuilder"`
- `filename`: HDF5 file containing the snapshot database
- `method`: Basis-construction method, currently `"POD"`

## Logged Output

During execution the solver reports:

- snapshot file name
- selected method
- number of snapshots
- state dimension
- number of resulting modes
- leading singular values
- captured energy for the first modes

## Singular Values and Energy

The POD basis is computed from the singular value decomposition of the snapshot
matrix:

```text
X^T = U S V^T
```

where:

- `U` contains the basis vectors written as `modes`
- `S` contains the singular values written as `eigenvals`

The reported energy is the cumulative fraction

```text
E(k) = (sigma_1^2 + ... + sigma_k^2) / (sigma_1^2 + ... + sigma_n^2)
```

This tells you how much of the snapshot content is represented by the first
`k` modes.

## Example

```text
solver =
{
  type = "ROMBasisBuilder";
  filename = "clamped.h5";
  method = "POD";
};
```

## Notes
- This solver does not collect snapshots itself; that is the job of `ROMSnapshotWriter`.
- In the ROM examples, `clamped_beam_prepare.pro` is the dedicated step that runs this solver after the training cases have populated the snapshot file.
- The resulting basis is consumed later by `ReducedOrderSolver`.