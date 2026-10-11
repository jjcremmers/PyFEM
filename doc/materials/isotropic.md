# Isotropic

`Isotropic` is a three-dimensional linear-elastic material model for
homogeneous isotropic solids. It can be used with small-strain or finite-strain
continuum elements, depending on the selected element formulation.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used. It is documentation support, not a substitute for verification of the
underlying equations or code.
```

## Overview

- **Material type:** `Isotropic`
- **Kinematics:** three-dimensional linear elasticity
- **Stress output:** `S11`, `S22`, `S33`, `S23`, `S13`, `S12`
- **Formulation:** total strain by default

The stress and strain vectors use the Voigt ordering

```text
[S11, S22, S33, S23, S13, S12]
```

The material computes stress from the strain vector and returns a constant
six-by-six elastic tangent matrix.

## Constitutive matrix

With

```{math}
\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)},
\qquad
\mu = \frac{E}{2(1+\nu)},
```

the constitutive relation is

```{math}
\boldsymbol{\sigma} = \boldsymbol{H}\,\boldsymbol{\varepsilon},
```

with

```{math}
\boldsymbol{H} =
\begin{bmatrix}
\lambda+2\mu & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2\mu & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2\mu & 0 & 0 & 0 \\
0 & 0 & 0 & \mu & 0 & 0 \\
0 & 0 & 0 & 0 & \mu & 0 \\
0 & 0 & 0 & 0 & 0 & \mu
\end{bmatrix}.
```

## Input parameters

```text
material =
{
  type = "Isotropic";
  E    = 100.0;
  nu   = 0.3;
};
```

| Parameter | Description | Type |
| --- | --- | --- |
| `type` | Must be `"Isotropic"`. | String |
| `E` | Young's modulus. | Float |
| `nu` | Poisson's ratio. | Float |

Use a consistent unit system for `E` and the resulting stresses.

## Examples

- [3D patch test](../../examples/ch02/PatchTest8_3D.pro)
- [3D Newton–Raphson example](../../examples/ch03/NewtonRaphson3D.pro)
- [3D contact example](../../examples/contact/contact_test3D_01.pro)
- [3D SLS cantilever](../../examples/elements/sls/sls_cantilever01.pro)
