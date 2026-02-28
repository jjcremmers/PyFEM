# DynEigSolver

`DynEigSolver` computes a selected number of eigenfrequencies and modes from `K` and `M` and reports the frequencies (rad/s and Hz). Eigenvectors are stored in the global data for postprocessing.

## Overview
- **Solver type:** `DynEigSolver`
- Problem: solves `K v = λ M v` and reports `ω = sqrt(λ)`
- Outputs: `eigenvecs` and `eigenvals` (rad/s)
- Reporting: logs mode number, angular frequency, and frequency in Hz

## Parameters
### Mandatory
- `type`: Must be set to `"DynEigSolver"`
