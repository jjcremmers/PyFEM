# RiksSolver

The `RiksSolver` (arc-length method) follows equilibrium paths through limit points by augmenting the Newton–Raphson iterations with an additional constraint on arc-length. It adapts step size to balance convergence effort.

## Overview
- **Solver type:** `RiksSolver`
- Method: arc-length with predictor–corrector iterations
- Control: maintains a constraint on combined displacement and load increments
- Step size: adapts `factor` per cycle toward an optimal iteration count
- Termination: stops when `lam > maxLam` or a cycle cap is reached

## Parameters
### Mandatory
- `type`: Must be set to `"RiksSolver"`
