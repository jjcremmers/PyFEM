# Modules

PyFEM provides solver and input/output modules that manage the execution of
analyses and the exchange of model and results data.

## Introduction

The modules described in this chapter are configured through PyFEM input files
and can be combined to build complete analyses. Solvers advance the solution,
while I/O modules read input data and write results for post-processing.

## Solvers

PyFEM provides nonlinear and other solver drivers to advance analyses.

```{toctree}
:maxdepth: -1

solvers/BuckEigSolver.md
solvers/DynEigSolver.md
solvers/DissipatedEnergySolver.md
solvers/ExplicitSolver.md
solvers/LinearSolver.md
solvers/ModalSolver.md
solvers/MultiSolver.md
solvers/NonlinearSolver.md
solvers/ROMBasisBuilder.md
solvers/ReducedOrderSolver.md
solvers/RiksSolver.md
solvers/StaggeredSolver.md
```

## I/O Modules

I/O modules read model data and write simulation results for post-processing.

```{toctree}
:maxdepth: -1

io/ContourWriter.md
io/DataDump.md
io/GraphWriter.md
io/HDF5Writer.md
io/InputReader.md
io/MeshWriter.md
io/OutputWriter.md
io/ROMSnapshotWriter.md
```
