# ModalSolver

The `ModalSolver` computes eigenmodes via the generalized eigenvalue problem for the assembled stiffness and mass matrices. It stores eigenvectors and values in the global data and terminates.

## Overview
- **Solver type:** `ModalSolver`
- Problem: solves `K v = λ M v`
- Outputs: eigenvectors and eigenvalues
- Termination: deactivates model after solve

## Parameters
### Mandatory
- `type`: Must be set to `"ModalSolver"`
