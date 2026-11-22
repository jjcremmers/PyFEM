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
solid bodies under small- and finite-strain assumptions. They are suitable for
bulk materials in two or three dimensions, including axisymmetric and
thermo-mechanical problems. The available formulations cover small- and
finite-strain kinematics, thermal expansion and conduction-related coupling,
and surface heat-transfer effects. Choose an element according to the assumed
strain regime, geometry, field variables, and material behavior. The individual
pages provide the supported degrees of freedom, integration details, required
material properties, and example `.pro` definitions.

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
They are intended for slender members and surface-like components where a
reduced structural model is more efficient than a full continuum mesh. PyFEM
includes linear and nonlinear beam elements, trusses, Kirchhoff and
Reissner–Mindlin plates and shells, and solid-like shell formulations. These
elements can represent axial deformation, bending, torsion, transverse shear,
and geometric nonlinearity. Select a formulation according to the component's
geometry, kinematic assumptions, and required degrees of freedom. Consult the
individual element pages for input properties and example `.pro` definitions.

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
These elements represent localized behavior that is not efficiently described
by standard continuum or structural elements. Interface and cohesive elements
can model separation, contact-like interactions, and damage along surfaces or
between material regions, while spring elements provide discrete stiffness or
connector behavior. They are often combined with other element types to model
joints, supports, delamination, or other localized mechanical effects. Consult
the individual pages for connectivity, degrees of freedom, and input syntax.

```{toctree}
:maxdepth: 1

interface.md
phasefield.md
spring.md
```
