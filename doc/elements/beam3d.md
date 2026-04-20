# Beam3D

The `Beam3D` element is a two-node, three-dimensional beam element for
large-displacement and large-rotation structural analysis.

## Overview

- **Element type:** `Beam3D`
- **Degrees of freedom per node:** `u`, `v`, `w`, `rx`, `ry`, `rz`
- **Supports:**
  - axial deformation
  - bending about two local axes
  - torsion
  - transverse shear deformation
  - large rotations with nonlinear solution procedures

`Beam3D` is currently implemented as a straight two-node beam with a frozen
reference geometry and a corotational-style local frame. The element supports
large rotations, but it should not be interpreted as a fully general
geometrically exact beam with arbitrary section coupling.

## Local Axes

The local beam `z` axis is aligned with the element axis from node 1 to node 2.
The local `x` and `y` axes define the section orientation.

The preferred way to prescribe the section orientation is through:

```text
orientation = [ x1, x2, x3 ];
```

This vector defines the preferred local `x` direction. It must:

- contain three components
- be non-zero
- not be parallel to the beam axis

If `orientation` is not given, `Beam3D` falls back to an internally constructed
seed direction:

- global `z` by default
- global `x` when the beam axis is nearly parallel to global `z`

## Material And Section Properties

### Standard scalar input

The legacy scalar input path uses:

```text
E, G, A, Ix, Iy, J, rho
```

These are interpreted as:

- `E`: Young's modulus
- `G`: shear modulus
- `A`: cross-sectional area
- `Ix`: second moment used for local `z`-direction bending stiffness
- `Iy`: second moment used for local `y`-direction bending stiffness
- `J`: torsional constant
- `rho`: mass density

### Explicit scalar stiffness input

The element also accepts direct section stiffness entries:

```text
EA, GAy, GAz, EIy, EIz, GJ
```

This is the recommended scalar interface when the section properties are already
available in stiffness form.

If `GAy` or `GAz` is not provided, the element uses:

- `ky * G * A` when `ky` is given
- `kz * G * A` when `kz` is given
- otherwise `5/6 * G * A`

### Section stiffness matrix

`Beam3D` also supports a direct section stiffness matrix:

```text
sectionStiffness = [ ... 36 values ... ];
```

The values are read as a row-major `6 x 6` matrix.

The internal ordering is:

```text
[ shear-x, shear-y, axial, bend-x, bend-y, torsion ]
```

Current limitation:

- coupling inside the first `3 x 3` block is supported
- coupling inside the second `3 x 3` block is supported
- coupling between these two blocks is not supported by the current
  reduced/full integration scheme

If off-diagonal coupling is present between the strain and curvature blocks, the
element raises a runtime error.

## Mass And Rotary Inertia

The legacy inertial input path uses:

```text
rho, A, Ix, Iy, J
```

from which the element constructs:

- mass per unit length: `rho * A`
- rotary inertia terms: `rho * Ix`, `rho * Iy`, `rho * J`

These values can also be supplied explicitly:

```text
massPerLength
rotaryInertia1
rotaryInertia2
polarRotaryInertia
```

The aliases below are accepted as well:

```text
rhoA, rhoIx, rhoIy, rhoJ
```

`massPerLength` must be positive. Rotary inertia values must be non-negative.

## Example

```text
BeamElem =
{
  type = "Beam3D";

  E           = 2.10e5;
  G           = 8.00e4;
  A           = 1.00e-2;
  Ix          = 1.00e-4;
  Iy          = 1.00e-4;
  J           = 2.00e-4;
  rho         = 7.85e-9;
  orientation = [ 0.0, 0.0, 1.0 ];
};
```

For direct section input:

```text
BeamElem =
{
  type = "Beam3D";

  orientation      = [ 0.0, 0.0, 1.0 ];
  sectionStiffness = [
    GAy, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, GAz, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, EA,  0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, EIy, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, EIz, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, GJ
  ];
};
```

## Examples

Beam3D example problems are available in:

```text
examples/elements/beam/beam3D/
```

The current example set includes:

- axial cantilever response
- transverse cantilever bending
- cantilever torsion
- vertical-member verification
- a deployable-ring demonstrator inspired by Goto et al. (1992)

## Current Scope And Limitations

- The element uses a two-node straight beam geometry.
- The reference geometry and reference triad are fixed after initialization.
- A prescribed section orientation is supported through `orientation`.
- A general uncoupled `6 x 6` section stiffness description is supported.
- Full coupling between translational strain resultants and curvature resultants
  is not supported.

For problems that depend strongly on general section coupling or more general
geometrically exact beam kinematics, this implementation should be treated as a
practical nonlinear 3D beam element rather than a complete general beam theory
framework.
