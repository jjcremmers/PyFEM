# Elements

PyFEM provides a comprehensive library of finite element formulations for
structural and solid mechanics analysis. Elements define the kinematic
assumptions, interpolation functions, and integration schemes that transform
continuum mechanics equations into discrete algebraic systems.

Elements are defined by creating named element groups in the `.pro` file. Each
group specifies the element type and its material properties.

See documentation for configuration examples and details.

## Continuum Elements

Continuum formulations describe the deformation and constitutive response of
solid bodies under small- and finite-strain assumptions.

```{toctree}
:maxdepth: 1

finitestrainaxisym.md
finitestraincontinuum.md
smallstrainaxisym.md
smallstraincontinuum.md
thermoaxiym.md
thermocontinuum.md
thermosmallstrainaxiym.md
thermosmallstraincontinuum.md
thermosurface.md
```

## Structural Elements

Structural formulations model beams, trusses, plates, and shell structures.
They describe slender and surface-like components efficiently, without
requiring a full continuum mesh. The library includes linear and nonlinear
beams, trusses, Kirchhoff and Reissner–Mindlin plates and shells, and solid-like
shell formulations. These
elements support common structural responses such as axial deformation,
bending, torsion, transverse shear, and geometric nonlinearity. The choice of
formulation depends on the intended kinematics, thickness-to-length ratio,
and whether transverse shear or finite rotations are important. See the
individual element pages for degrees of freedom, required input properties,
integration assumptions, and example `.pro` definitions.

```{toctree}
:maxdepth: 1

beam3d.md
beamlin.md
beamnl.md
kirchhoffbeam.md
plate.md
reissnermindlinshell.md
sls.md
solidlikeshell.md
timoshenkobeam.md
truss.md
```

## Special Elements

Special formulations include interface, cohesive-zone, and spring elements.

```{toctree}
:maxdepth: -1

interface.md
phasefield.md
spring.md
```
