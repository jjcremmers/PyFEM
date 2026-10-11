# PlaneStrain

`PlaneStrain` is a linear-elastic material model for two-dimensional plane
strain analyses. It is appropriate for thick bodies or long structures where
the out-of-plane strain can be assumed negligible.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used. It is documentation support, not a substitute for verification of the
underlying equations or code.
```

## Overview

- **Material type:** `PlaneStrain`
- **Kinematics:** small-strain, two-dimensional plane strain
- **Stress output:** `S11`, `S22`, `S12`

Under the plane strain assumption, the out-of-plane strain is constrained to
zero. The material returns the in-plane stress vector

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
\dfrac{E(1-\nu)}{(1+\nu)(1-2\nu)} &
\dfrac{E\nu}{(1+\nu)(1-2\nu)} & 0 \\
\dfrac{E\nu}{(1+\nu)(1-2\nu)} &
\dfrac{E(1-\nu)}{(1+\nu)(1-2\nu)} & 0 \\
0 & 0 & \dfrac{E}{2(1+\nu)}
\end{bmatrix}.
```

The model returns this constant matrix as its tangent.

## Input parameters

The material is defined inside an element material block:

```text
material =
{
  type = "PlaneStrain";
  E    = 100.0;
  nu   = 0.3;
};
```

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"PlaneStrain"`. | String | Required |
| `E` | Young's modulus. | Float | Required |
| `nu` | Poisson's ratio. | Float | Required |

Use a consistent unit system for `E` and the resulting stresses. The selected
element must provide a compatible two-dimensional strain vector.
