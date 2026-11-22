# MultiSolver

`MultiSolver` runs a sequence of solver blocks in order. When one solver finishes (deactivates the model), `MultiSolver` activates the next solver and continues until all are complete.

## Overview
- **Solver type:** `MultiSolver`
- Orchestration: imports and constructs solver classes by name
- Progression: advances to the next solver when `globdat.active` becomes false

## Parameters
### Mandatory
- `type`: Must be set to `"MultiSolver"`
- `solvers`: List of solver block names to run in order

For each name in `solvers`, define a block with at least a `type`.
