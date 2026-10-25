# Reissner–Mindlin shell

## Background

The `ReissnerMindlinShell` element models flat and curved shell midsurfaces
with transverse shear deformation and geometrically nonlinear kinematics. It
is a four-node quadrilateral element with six degrees of freedom per node:

```text
u, v, w, rx, ry, rz
```

The translational degrees of freedom are `u`, `v`, and `w`; `rx` and `ry` are
the director rotations, and `rz` is the drilling rotation. The element is
formulated in a local shell frame and supports layer-wise integration through
the thickness.

:::{warning}
Only four-node quadrilateral midsurfaces are supported. The element can be
used for flat and initially curved geometries, but other interpolation types
are not supported.
:::

## Implementation

The implementation is in `pyfem/elements/ReissnerMindlinShell.py`. It uses
the same material and laminate definitions as the [`Plate`](plate.md)
element.

At each integration point, the element evaluates the current shell basis and
director, then computes the generalized membrane, bending, and transverse
shear strains. The material response is evaluated layer by layer. In compact
form, the in-plane resultants and bending moments follow the laminate relation

$$
\begin{bmatrix}
\boldsymbol{N} \\
\boldsymbol{M}
\end{bmatrix}
=
\begin{bmatrix}
\boldsymbol{A} & \boldsymbol{B} \\
\boldsymbol{B} & \boldsymbol{D}
\end{bmatrix}
\begin{bmatrix}
\boldsymbol{\varepsilon}^{0} \\
\boldsymbol{\kappa}
\end{bmatrix},
$$

while transverse shear uses the layer shear stiffness and the shear-correction
factor. The tangent stiffness contains both material and geometric
contributions. The director Jacobian is evaluated locally, and the tangent is
obtained using the perturbation specified by `tangentPerturbation`.

By default, the transverse shear terms use selective reduced integration. A
small drilling stabilization prevents the `rz` degree of freedom from being
singular. The default drilling scale is deliberately small so that it does
not significantly affect the shell response.

Stress output is reported at the top and bottom surfaces using these labels:

```text
S11_top, S22_top, S12_top, S13_top, S23_top,
S11_bot, S22_bot, S12_bot, S13_bot, S23_bot
```

## Input parameters

The element block must contain `type = "ReissnerMindlinShell"` and either a
single material definition or a laminate definition.

### Single material

This is the form used by
[clamped_bending.pro](../../examples/elements/reissnermindlin/linear/clamped_bending.pro):

```text
ShellElem =
{
  type = "ReissnerMindlinShell";

  material =
  {
    E   = 2.1e5;
    nu  = 0.3;
    rho = 7.85e-9;
  };

  thickness    = 1.0;
  drillingScale = 1.0e-6;
};
```

| Parameter | Description | Optional | Default value |
| --- | --- | --- | --- |
| `type` | Must be `"ReissnerMindlinShell"`. | No | — |
| `material` | Single material definition. | No | — |
| `thickness` | Shell thickness. | No | — |
| `E` | Young's modulus; a scalar gives equal principal moduli. | No | — |
| `nu` or `nu12` | Poisson's ratio. | No | — |
| `rho` | Density, used for inertia. | No | Required |
| `G12` | In-plane shear modulus. | Yes | $E/[2(1+\nu)]$ |
| `G13` | 1–3 shear modulus. | Yes | `G12` |
| `G23` | 2–3 shear modulus. | Yes | `G12` |
| `shearCorrection` | Transverse shear correction factor. | Yes | $5/6$ |
| `tangentPerturbation` | Perturbation used to compute the numerical tangent. | Yes | $10^{-7}$ |
| `drillingScale` | Stabilization scale for `rz`. | Yes | $10^{-6}$ |
| `reducedShearIntegration` | Enables selective reduced integration of shear terms. | Yes | `true` |

### Laminate

For a layered shell, define the material names in `materials` and the layer
names in `layers`. This is demonstrated by
[curved_cantilever_composite.pro](../../examples/elements/reissnermindlin/curved_cantilever_composite.pro).

```text
ShellElem =
{
  type      = "ReissnerMindlinShell";
  materials = ["UD"];
  layers    = ["ply0_bot", "ply90_bot", "ply90_top", "ply0_top"];

  UD =
  {
    E1   = 1.35e5;
    E2   = 1.0e4;
    nu12 = 0.3;
    G12  = 5.0e3;
    G13  = 4.0e3;
    G23  = 3.8e3;
    rho  = 1.6e-9;
  };

  ply0_bot  = { material = "UD"; theta = 0.0;  thickness = 0.005; };
  ply90_bot = { material = "UD"; theta = 90.0; thickness = 0.005; };
  ply90_top = { material = "UD"; theta = 90.0; thickness = 0.005; };
  ply0_top  = { material = "UD"; theta = 0.0;  thickness = 0.005; };
};
```

The element-level parameters are:

| Parameter | Description | Optional | Default value |
| --- | --- | --- | --- |
| `type` | Must be `"ReissnerMindlinShell"`. | No | — |
| `materials` | Names of material blocks. | No | — |
| `layers` | Ordered layer-block names, from bottom to top. | No | — |
| `shearCorrection` | Transverse shear correction factor. | Yes | $5/6$ |
| `tangentPerturbation` | Perturbation used to compute the numerical tangent. | Yes | $10^{-7}$ |
| `drillingScale` | Stabilization scale for `rz`. | Yes | $10^{-6}$ |
| `reducedShearIntegration` | Enables selective reduced integration of shear terms. | Yes | `true` |

Material and layer blocks use the same parameters as the plate laminate:

| Parameter | Description | Optional | Default value |
| --- | --- | --- | --- |
| `E1` and `E2`, or `E` | Principal Young's moduli; `E` may be scalar or a two-entry list. | No | — |
| `nu12` or `nu` | Poisson's ratio in material axes. | No | — |
| `rho` | Material density. | No | Required |
| `G12` | In-plane shear modulus. | Yes | $E_1/[2(1+\nu_{12})]$ |
| `G13` | 1–3 shear modulus. | Yes | `G12` |
| `G23` | 2–3 shear modulus. | Yes | `G12` |
| `material` | Name of a material in `materials`. | No | — |
| `thickness` | Thickness of the layer. | No | — |
| `theta` | Layer orientation in degrees. | Yes | $0^\circ$ |

The total shell thickness is the sum of the layer thicknesses:

$$
t = \sum_{i=1}^{n} t_i.
$$

## Examples

The examples in `examples/elements/reissnermindlin` cover the main use cases:

- [Clamped bending](../../examples/elements/reissnermindlin/linear/clamped_bending.pro)
  demonstrates a linear analysis of a single-material shell.
- [Curved composite cantilever](../../examples/elements/reissnermindlin/curved_cantilever_composite.pro)
  demonstrates a curved shell with four orthotropic layers.
- [Pinched hemisphere](../../examples/elements/reissnermindlin/nonlinear/pinched_hemisphere.pro)
  demonstrates geometrically nonlinear analysis.
- [Panel eigenfrequencies](../../examples/elements/reissnermindlin/dynamic/panel_eigenfrequencies.pro)
  demonstrates dynamic analysis using shell inertia.

The shell can be used with a linear, nonlinear, or eigenvalue solver. The
appropriate solver and output modules are selected in each example.
