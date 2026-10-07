# Plate

`Plate` is a flat, small-displacement plate element for meshes in the global
`x-y` plane. Each node has five degrees of freedom: `u`, `v`, `w`, `rx`, and
`ry`.

:::{warning}
The plate element supports 4-node elements. Other surface elements (3-, 6-,
and 8-node elements) can be used subject to the shape functions used by the
mesh, but these configurations have not been tested.
:::

The element accepts either a single material or a laminate. Layer angles are
given in degrees. Layers are integrated about the laminate mid-plane, so their
ordered thicknesses define the total plate thickness.

## Implementation

The implementation is in `pyfem/elements/Plate.py` and uses `Laminate` from
`pyfem/elements/Composite.py`. At each integration point it computes membrane
strains from `u` and `v`, curvatures from `rx` and `ry`, and transverse shear
strains from `w`, `rx`, and `ry`.

The membrane resultants and bending moments use the laminate constitutive
matrices:

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
\end{bmatrix}.
$$

`A` is the extensional stiffness, `B` is the membrane-bending coupling
stiffness, and `D` is the bending stiffness. Transverse shear uses the
laminate shear matrix and a default shear-correction factor of
$k_s = \frac{5}{6}$.

The element also provides a consistent mass matrix. Mass and rotary inertia
are calculated from the layer densities and layer positions; therefore `rho`
is required for dynamic analyses.

Stress output is evaluated at the bottom and top surfaces with these labels:

```text
s11bot, s22bot, s12bot, s11top, s22top, s12top
```

## Input parameters

The element block must contain `type = "Plate"` and one of the following
material definitions.

### Single isotropic layer: parameter overview

This is the compact form used by
[plate_cantilever01.pro](../../examples/elements/plate/plate_cantilever01.pro).

| Parameter | Location | Required? | Default | Description |
| --- | --- | --- | --- | --- |
| `type` | element block | Mandatory | — | Must be `"Plate"`. |
| `material` | element block | Mandatory | — | Single material definition. |
| `E` | `material` | Mandatory | — | Young's modulus; a scalar gives equal principal moduli. |
| `nu` or `nu12` | `material` | Mandatory | — | Poisson's ratio. |
| `rho` | `material` | Mandatory | — | Density, used for the mass matrix. |
| `thickness` | element block | Mandatory | — | Plate thickness. |
| `G12` | `material` | Optional | $E/[2(1+\nu)]$ | In-plane shear modulus. |
| `G13` | `material` | Optional | `G12` | 1–3 shear modulus. |
| `G23` | `material` | Optional | `G12` | 2–3 shear modulus. |
| `theta` | element block | Optional | $0^\circ$ | Material orientation in degrees. |
| `shearCorrection` | element block | Optional | $5/6$ | Transverse shear correction factor. |

### Single isotropic material

Use `material` and `thickness` for a single-layer plate:

```text
PlateElem =
{
  type = "Plate";

  material =
  {
    E       = 1.0e6;
    nu      = 0.25;
    rho     = 1.0e3;
  };

  thickness = 0.1;
};
```

The material properties consumed by the plate are:

- `E`: Young's modulus.
- `nu` `: Poisson's ratio .
- `rho`: density.
- `thickness`: plate thickness.

For this form, the layer angle defaults to `0` degrees.

### Laminate

For multiple materials or layers, define the material names in `materials`
and the layer names in `layers`:

### Multiple layers: parameter overview

This is the form used by
[plate_cantilever02.pro](../../examples/elements/plate/plate_cantilever02.pro).
The `layers` list defines the order through the thickness.

| Parameter | Location | Required? | Default | Description |
| --- | --- | --- | --- | --- |
| `type` | element block | Mandatory | — | Must be `"Plate"`. |
| `materials` | element block | Mandatory | — | Names of material blocks. |
| `layers` | element block | Mandatory | — | Ordered layer-block names, from bottom to top. |
| `E1` and `E2`, or `E` | material block | Mandatory | — | Principal Young's moduli; `E` may be scalar or a two-entry list. |
| `nu12` or `nu` | material block | Mandatory | — | Poisson's ratio in material axes. |
| `rho` | material block | Mandatory | — | Material density. |
| `G12` | material block | Optional | $E_1/[2(1+\nu_{12})]$ | In-plane shear modulus. |
| `G13` | material block | Optional | `G12` | 1–3 shear modulus. |
| `G23` | material block | Optional | `G12` | 2–3 shear modulus. |
| `material` | layer block | Mandatory | — | Name of a material in `materials`. |
| `thickness` | layer block | Mandatory | — | Thickness of the layer. |
| `theta` | layer block | Optional | $0^\circ$ | Layer orientation in degrees. |
| `shearCorrection` | element block | Optional | $5/6$ | Transverse shear correction factor. |

The total laminate thickness is

$$
t = \sum_{i=1}^{n} t_i.
$$

```text
PlateElem =
{
  type      = "Plate";
  materials = ["UD"];
  layers    = ["l0", "l90", "l0"];

  UD =
  {
    E1   = 1.0e6;
    E2   = 5.0e5;
    nu12 = 0.25;
    G12  = 4.0e5;
    rho  = 1.0e3;
  };

  l0 = { material = "UD"; theta = 0.0;  thickness = 0.05; };
  l90 = { material = "UD"; theta = 90.0; thickness = 0.05; };
};
```

Each layer requires `material` and `thickness`; `theta` is optional and
defaults to `0` degrees. The order in `layers` is the order through the
thickness, from bottom to top. Total thickness is the sum of the layer
thicknesses:

$$
t = \sum_{i=1}^{n} t_i.
$$

The optional element property `shearCorrection` overrides the default `5/6`:

```text
shearCorrection = 0.8333333333;
```

Some older examples contain a `stack` entry. The current implementation does
not read `stack`; use the ordered `layers` list to define the laminate.

## Examples

The two files in `examples/elements/plate` are the reference examples for this
element:

- [Single-material cantilever](../../examples/elements/plate/plate_cantilever01.pro)
  uses `material` with isotropic properties (`E`, `nu`, and `rho`) and a
  plate-level `thickness` of `0.1`.
- [Layered cantilever](../../examples/elements/plate/plate_cantilever02.pro)
  uses one orthotropic material (`UD`) and the layer sequence
  `l0`–`l90`–`l0`, with each layer assigned its own angle and thickness.

Both examples use a `LinearSolver`, write VTK output with `MeshWriter`, and
enable screen output with `OutputWriter`. The essential single-material input
from the first example is:

```text
input = "plate_cantilever01.dat";

PlateElem =
{
  type = "Plate";
  material = { E = 1.0e6; nu = 0.0; rho = 1.0e3; };
  thickness = 0.1;
};

solver = { type = "LinearSolver"; };
```
