# Elements Overview

PyFEM provides a comprehensive library of finite element formulations for structural and solid mechanics analysis. Elements define the kinematic assumptions, interpolation functions, and integration schemes that transform continuum mechanics equations into discrete algebraic systems.

## Element Types
- **Continuum elements:** 2D, 3D, and axisymmetric solids
- **Structural elements:** Beams, trusses, plates, shells
- **Interface elements:** Cohesive zone models
- **Special elements:** Springs and connectors

## Configuration
Elements are defined by creating named element groups in the `.pro` file. Each group specifies the element type and its material properties.

## Available Element Pages

```{toctree}
:maxdepth: 1

beamnl.md
finitestrainaxisym.md
finitestraincontinuum.md
interface.md
kirchhoffbeam.md
plate.md
sls.md
smallstrainaxisym.md
smallstraincontinuum.md
spring.md
timoshenkobeam.md
truss.md
```

See documentation for configuration examples and details.
