# Reduced-Order Modeling Examples

This directory contains a compact reduced-order modeling workflow for a clamped
beam problem. The files cover both the full-order model (FOM) runs used to
generate snapshot data and the reduced-order model (ROM) runs that reuse a
precomputed modal basis.

## Use Cases

| Use case | Recommended input file | Purpose |
| --- | --- | --- |
| Run a reference full-order nonlinear simulation | `clamped_beam_FOM.pro` | Baseline solution on the original mesh and DOF space |
| Compute eigenmodes for a clamped beam | `clamped_beam.pro` | Uses `DynEigSolver` to generate modal information |
| Build the POD basis after snapshot collection | `clamped_beam_prepare.pro` | Reads the snapshot HDF5 file and writes the reduced basis |
| Collect training snapshots for offline ROM construction | `clamped_beam_training1.pro`, `clamped_beam_training2.pro`, `clamped_beam_training3.pro` | Runs several loading cases and appends state snapshots to the ROM file |
| Run a reduced model from an existing mode file | `clamped_beam_ROM.pro` | Uses `ReducedOrderSolver` with an HDF5 file containing a `modes` dataset |
| Run the alternate online ROM case | `clamped_beam_online.pro` | Same online idea, but with a different `modeCount` |

## Recommended Workflow

Use the files in this order when building a ROM from scratch:

1. Run `pyfem clamped_beam_training1.pro`
2. Run `pyfem clamped_beam_training2.pro`
3. Run `pyfem clamped_beam_training3.pro`
4. Run `pyfem clamped_beam_prepare.pro`
5. Run `pyfem clamped_beam_ROM.pro`

The training runs append snapshot states through `ROMSnapshotWriter`,
`clamped_beam_prepare.pro` then builds the reduced basis from that collected
snapshot database, and the ROM run uses the stored POD modes for
Newton iterations in reduced coordinates.

If you only want a baseline comparison, run `pyfem clamped_beam_FOM.pro`.

## File Roles

### Full-order runs

- `clamped_beam_FOM.pro`: nonlinear full-order reference problem on `clamped_beam.dat`
- `clamped_beam.pro`: eigenvalue-based beam analysis using `DynEigSolver`

### Offline ROM generation

- `clamped_beam_training1.pro`: training case based on `clamped_beam1.dat`
- `clamped_beam_training2.pro`: training case based on `clamped_beam2.dat`
- `clamped_beam_training3.pro`: training case based on `clamped_beam3.dat`
- `clamped_beam_prepare.pro`: pure basis-construction step using `ROMBasisBuilder` with `clamped.h5`

### Online ROM execution

- `clamped_beam_ROM.pro`: online ROM run using `ReducedOrderSolver` and `clamped.h5`
- `clamped_beam_online.pro`: online ROM run using `ReducedOrderSolver` and `clamped.h5`

## Important Note About the Prepare Step

`clamped_beam_prepare.pro` does not collect snapshots itself.

Its only job is to read the `state` snapshots already stored in `clamped.h5`
and write the POD basis (`modes`) and singular values (`eigenvals`) back to
that same file.

So the offline sequence is always:

1. collect snapshots with the training cases
2. build the basis with `clamped_beam_prepare.pro`
3. run the reduced model with `clamped_beam_ROM.pro` or `clamped_beam_online.pro`

## Data Files

- `clamped_beam.dat`: baseline beam input used by the FOM and one ROM run
- `clamped_beam1.dat`, `clamped_beam2.dat`, `clamped_beam3.dat`: training variants used to enrich the snapshot set
- `clamped.h5`: HDF5 snapshot and mode database produced by the offline workflow

## Expected Outputs

- `*.pvd` and `*.vtu`: visualization files written by `MeshWriter`
- `*.h5`: snapshot and mode storage used by the ROM pipeline

The ROM solver expects the HDF5 file to contain a dataset named `modes`. The
offline preparation step constructs that dataset using POD from the collected
state snapshots.

## Understanding the Basis-Builder Output

When `ROMBasisBuilder` runs, it reports singular values and captured energy.

- `leading sigma` and `second sigma` are the largest singular values of the snapshot matrix.
- Larger singular values indicate more important deformation patterns in the training data.
- `energy(1 mode)` is the fraction of total snapshot content captured by the first mode.
- `energy(2 modes)` is the fraction captured by the first two modes together.

In practice, these values help decide how many modes you should keep in
`ReducedOrderSolver`.

For example:

- if `energy(1 mode)` is already close to `1.0`, one mode may be enough for a useful ROM
- if the cumulative energy rises slowly, you should keep more modes
- if the first few singular values drop very quickly, the problem is strongly compressible

The energy values are computed from the squared singular values:

```text
E(k) = (sigma_1^2 + ... + sigma_k^2) / (sigma_1^2 + ... + sigma_n^2)
```

So they are best interpreted as cumulative variance or cumulative information
content of the snapshot set.

## Practical Guidance

- Use `clamped_beam_FOM.pro` when you need a reference response to compare with the reduced solution.
- Use the three training cases when one load case is not rich enough to build a robust reduced basis.
- Start with a small `modeCount` in `clamped_beam_ROM.pro` and increase it if the reduced solution is too inaccurate.
- If a ROM run fails immediately, check that the mode file exists and that its number of rows matches the number of global DOFs.