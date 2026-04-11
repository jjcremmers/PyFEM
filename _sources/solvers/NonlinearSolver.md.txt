# NonlinearSolver

The `NonlinearSolver` advances the analysis in cycles using a Newton–Raphson procedure. At each cycle it assembles the tangent stiffness and internal force, solves for the displacement increment, updates the state, and checks convergence against a residual norm. Loads can be defined via a time-based function or a tabulated sequence and may include multiple cases.

## Overview
- **Solver type:** `NonlinearSolver`
- Method: standard Newton–Raphson with residual norm check
- Termination: when convergence is reached each cycle; run stops at `maxCycle` or when `lam > maxLam`
- Load control: `loadFunc` evaluated over time or `loadTable` per cycle; additional `loadCases` supported

## Parameters
### Mandatory
- `type`: Must be set to `"NonlinearSolver"`
