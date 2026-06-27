# ExplicitSolver

The `ExplicitSolver` advances dynamics via a central-difference explicit time-integration scheme using a lumped mass matrix. It updates velocity and displacement, computes accelerations from the residual, and reports kinetic energy.

## Overview
- **Solver type:** `ExplicitSolver`
- Method: explicit central differences with half-step velocity updates
- Mass: uses lumped mass matrix from assembly
- Loads: scalar load factor `lam(t)` applied to `fhat`
- Termination: run stops at `maxCycle`

## Parameters
### Mandatory
- `type`: Must be set to `"ExplicitSolver"`
