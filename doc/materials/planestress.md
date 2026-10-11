# PlaneStress

`PlaneStress` is a linear-elastic material model for two-dimensional plane
stress analyses. It is appropriate for thin bodies where the out-of-plane
stress can be assumed negligible.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used. It is documentation support, not a substitute for verification of the
underlying equations or code.
```

## Overview

- **Material type:** `PlaneStress`
- **Kinematics:** small-strain, two-dimensional plane stress
- **Stress output:** `S11`, `S22`, `S12`

The material returns the in-plane stress vector

```{math}
\boldsymbol{\sigma} = [S_{11}, S_{22}, S_{12}]^T
```

from the engineering strain vector

```{math}
\boldsymbol{\varepsilon} = [\varepsilon_{11}, \varepsilon_{22}, \gamma_{12}]^T.
```

## Constitutive matrix

The stress is computed as

```{math}
\boldsymbol{\sigma} = \boldsymbol{H}\,\boldsymbol{\varepsilon},
```

where the constitutive matrix is

```{math}
\boldsymbol{H} =
\begin{bmatrix}
\dfrac{E}{1-\nu^2} & \dfrac{E\nu}{1-\nu^2} & 0 \\
\dfrac{E\nu}{1-\nu^2} & \dfrac{E}{1-\nu^2} & 0 \\
0 & 0 & \dfrac{E}{2(1+\nu)}
\end{bmatrix}.
```

The model returns this constant matrix as its tangent.

## Input parameters

The material is defined inside an element material block:

```text
material =
{
  type = "PlaneStress";
  E    = 100.0;
  nu   = 0.3;
};
```

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"PlaneStress"`. | String | |
| `E` | Young's modulus. | Float | |
| `nu` | Poisson's ratio. | Float | |

Use a consistent unit system for `E` and the resulting stresses. The selected
element must provide a compatible two-dimensional strain vector.

## Examples

The repository contains several examples that use `PlaneStress`:

- [Patch test 3](../../examples/ch02/PatchTest3.pro)
- [Patch test 8](../../examples/ch02/PatchTest8.pro)
- [Phase-field patch test](../../examples/elements/phasefield/PhaseFieldPatch4.pro)
- [Traction oscillations example](../../examples/ch13/TractionOscillations.pro)

These examples demonstrate the material in patch tests, continuum analyses,
and a phase-field application.
