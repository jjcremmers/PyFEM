# BuckEigSolver

The `BuckEigSolver` performs a buckling eigenvalue analysis. It solves for the prebuckling static state, assembles the updated stiffness, and computes buckling eigenvalues/eigenvectors.

## Overview
- **Solver type:** `BuckEigSolver`
- Prebuckling: solves `K0 a = f_ext` for the static state
- Buckling: computes eigenpairs from the generalized problem using the initial and updated stiffness
- Outputs: `eigenvals` (critical loads) and `eigenvecs` (buckling modes)

## Parameters
### Mandatory
- `type`: Must be set to `"BuckEigSolver"`
