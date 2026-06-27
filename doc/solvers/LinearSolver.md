# LinearSolver

The `LinearSolver` performs a single linear solve for the global displacement vector given the assembled tangent stiffness and external force. It then computes derived quantities and terminates the analysis.

## Overview
- **Solver type:** `LinearSolver`
- Method: single-step linear solution `K a = f_ext`
- Outputs: updates state, computes `Dstate`, internal forces, commits history
- Termination: deactivates model after the single solve

## Parameters
### Mandatory
- `type`: Must be set to `"LinearSolver"`
