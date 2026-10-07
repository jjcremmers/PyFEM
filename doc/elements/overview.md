# Elements Overview

PyFEM provides a comprehensive library of finite element formulations for structural and solid mechanics analysis. Elements define the kinematic assumptions, interpolation functions, and integration schemes that transform continuum mechanics equations into discrete algebraic systems.

## Configuration
Elements are defined by creating named element groups in the `.pro` file. Each group specifies the element type and its material properties.

## Available Element Models

```{toctree}
:maxdepth: 1

beam3d.md
beamnl.md
finitestrainaxisym.md
finitestraincontinuum.md
interface.md
kirchhoffbeam.md
plate.md
reissnermindlinshell.md
sls.md
smallstrainaxisym.md
smallstraincontinuum.md
spring.md
timoshenkobeam.md
truss.md
```

See documentation for configuration examples and details.
