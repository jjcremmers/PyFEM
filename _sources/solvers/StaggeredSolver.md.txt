# StaggeredSolver

`StaggeredSolver` splits the solution into two subproblems handled by configured sub-solvers with disjoint dof sets. It alternates solves for each subset within a cycle and can perform Newton iterations for sub-solvers of type `Nonlinear`.

## Overview
- **Solver type:** `StaggeredSolver`
- Sub-solvers: two configured blocks (`solver1`, `solver2`) with names, types, and `dofTypes` defining their scope
- Load control: shared `loadFunc` and optional `loadCases` applied via sub-solver constrainers
- Termination: when reaching `maxCycle`

## Parameters
### Mandatory
- `type`: Must be set to `"StaggeredSolver"`
- `solver1` / `solver2`: Sub-blocks describing each solver
