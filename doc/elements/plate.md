# Plate

The `Plate` element models flat plates in the global `x-y` plane under small
displacements. Each node has five degrees of freedom: `u`, `v`, `w`, `rx`,
and `ry`.

:::{admonition} Element support
:class: plate-support-warning
The plate element supports 4-node elements. Other surface elements (3-, 6-,
and 8-node elements) can be used subject to the shape functions used by the
mesh, but these configurations have not been tested.
:::

The element supports both a single material layer and a multilayer laminate.
Layer angles are specified in degrees. Layers are integrated about the
laminate mid-plane; their ordered thicknesses determine the total plate
thickness.

## Implementation

The implementation is in `pyfem/elements/Plate.py` and uses the `Laminate`
class from `pyfem/elements/Composite.py`. At each integration point, it
computes membrane strains from `u` and `v`, curvatures from `rx` and `ry`, and
transverse shear strains from `w`, `rx`, and `ry`.

The membrane resultants and bending moments use the laminate constitutive
matrices:

```{math}
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
```

Here, `A` is the extensional stiffness, `B` is the membrane-bending coupling
stiffness, and `D` is the bending stiffness. Transverse shear uses the
laminate shear matrix with the default correction factor
$k_s = \frac{5}{6}$.

The element also provides a consistent mass matrix. Mass and rotary inertia
are calculated from the layer densities and layer positions, so `rho` must be
provided when dynamic effects are included.

Stress output is evaluated at the bottom and top surfaces with these labels:

```text
s11bot, s22bot, s12bot, s11top, s22top, s12top
```

## Input parameters

The element block must contain `type = "Plate"` and one of the following
material definitions.

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

The element-level parameters for a single isotropic layer are:

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type`      | Must be `"Plate"`. | String | — |
| `material`  | Single material definition. | Data block | — |
| `thickness` | Plate thickness. | Float    | — |
| `shearCorrection` | Transverse shear correction factor. | Float | Optional, default = $5/6$ |

The material definition is:

| Parameter | Description | Type |  Remarks  |
| ---  | --- | --- | --- |
| `E`  | Young's modulus; a scalar gives equal principal moduli. | Float | — |
| `nu` | Poisson's ratio. | Float | — |
| `rho` | Density, used for the mass matrix. | Float | Required |


### Laminate

For multiple materials or layers, define the material names in `materials`
and the layer names in `layers`:

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

  l0  = { material = "UD"; theta = 0.0;  thickness = 0.05; };
  l90 = { material = "UD"; theta = 90.0; thickness = 0.05; };
};
```

The element-level parameters are:

| Parameter | Description | Type |  Remarks  |
| ---  | --- | --- | --- |
| `type` | Must be `"Plate"`. | String | — |
| `materials` | Names of material blocks. | List of strings | — |
| `layers` | Ordered layer-block names, from bottom to top. | List of strings | — |
| `shearCorrection` | Transverse shear correction factor. | Float | Optional, default = $5/6$ |

Define one material block for each name in the `materials` list:

| Parameter | Description | Type |  Remarks  |
| ---  | --- | --- | --- |
| `E1` and `E2`, or `E` | Principal Young's moduli; `E` may be scalar or a two-entry list. | Float | — |
| `nu12` or `nu` | Poisson's ratio in material axes. | Float | — |
| `rho` | Material density. | Float | Optional, default = 0.0 |
| `G12` | In-plane shear modulus | Float | Optional, default = $E_1/[2(1+\nu_{12})]$ |
| `G13` | 1–3 shear modulus. | Float | Optional, default = `G12` |
| `G23` | 2–3 shear modulus. | Float | Optional, default = `G12` |

Define one layer block for each name in the `layers` list:

| Parameter | Description | Type |  Remarks  |
| ---  | --- | --- | --- |
| `material` | Name of a material in `materials`. | String | — |
| `thickness` | Thickness of the layer. | Float | — |
| `theta` | Layer orientation in degrees. | Float | Optional, default = $0^\circ$ |


The order in `layers` is the order through the thickness, from bottom to top.
The total laminate thickness is the sum of the individual layer thicknesses:

```{math}
t = \sum_{i=1}^{n} t_i.
```

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
