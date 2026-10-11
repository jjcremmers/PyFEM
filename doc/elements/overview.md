# Elements

## Introduction

PyFEM provides a comprehensive library of finite element formulations for
structural and solid mechanics analysis. Elements define the kinematic
assumptions, interpolation functions, and integration schemes that transform
continuum mechanics equations into discrete algebraic systems.

Elements are defined by creating named element groups in the `.pro` file. Each group specifies the element type and its material properties.

See documentation for configuration examples and details.

## Continuum Elements

Continuum formulations describe the deformation and constitutive response of
solid bodies under small- and finite-strain assumptions.

```{toctree}
:maxdepth: -1

finitestrainaxisym.md
finitestraincontinuum.md
smallstrainaxisym.md
smallstraincontinuum.md
```

## Structural Elements

Structural formulations model beams, trusses, plates, and shell structures.

```{toctree}
:maxdepth: -1

beam3d.md
beamnl.md
kirchhoffbeam.md
plate.md
reissnermindlinshell.md
sls.md
timoshenkobeam.md
truss.md
```

## Special Elements

Special formulations include interface, cohesive-zone, and spring elements.

```{toctree}
:maxdepth: -1

interface.md
spring.md
```
