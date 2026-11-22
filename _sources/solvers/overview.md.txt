# Solvers

Solvers control how PyFEM advances a finite element analysis from one state to
the next. The solver framework includes linear and nonlinear static methods,
explicit and dynamic analyses, eigenvalue and modal calculations, arc-length
continuation, staggered multiphysics schemes, and reduced-order workflows.

Select a solver
according to the governing equations, physical problem, and desired solution
strategy. Solver settings are defined in the `.pro` input file. The pages below
describe the available parameters, supported workflows, use cases, and
examples for configuring each solver.

```{toctree}
:maxdepth: 1

BuckEigSolver.md
DynEigSolver.md
DissipatedEnergySolver.md
ExplicitSolver.md
LinearSolver.md
ModalSolver.md
MultiSolver.md
NonlinearSolver.md
ROMBasisBuilder.md
ReducedOrderSolver.md
RiksSolver.md
StaggeredSolver.md
```
