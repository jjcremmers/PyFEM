# SLS (solid-like shell)

## Background

The `SLS` element is a solid-like shell formulation for layered and laminated
structures. It is intended for composite plates, sandwich panels, and other
shells for which the through-thickness material response is important.

Unlike conventional shell elements, the SLS formulation uses only the three
translational degrees of freedom at the external nodes:

```text
u, v, w
```

Additional kinematic variables are introduced internally and statically
condensed at element level. The resulting global system therefore contains
only the translational degrees of freedom, while still representing bending,
transverse shear, and through-thickness behavior.

:::{warning}
The current implementation supports only 8-node and 16-node SLS elements.
The 8-node formulation uses an assumed-natural-strain treatment; the 16-node
formulation uses the higher-order interpolation without that treatment.
:::

`SLS` is a backward-compatible name for the `SolidLikeShell` implementation.
It supports both geometrically nonlinear analyses and material models with
incremental behavior.

## Implementation

The compatibility class is defined in `pyfem/elements/SLS.py`; the element
implementation is in `pyfem/elements/SolidLikeShell.py`. The element performs
the following operations at each integration point:

1. Construct the shell kinematics from the current nodal translations.
2. Evaluate the material response at each through-thickness layer point.
3. Assemble material and geometric stiffness contributions.
4. Condense the internal element degrees of freedom before returning the
   element stiffness and internal force.

For a layer at thickness coordinate $z$, the generalized shell strain is
represented in the usual membrane-bending form:

```{math}
\boldsymbol{\varepsilon}(z)
=
\boldsymbol{\varepsilon}^{0}
+ z\,\boldsymbol{\kappa},
```

where $\boldsymbol{\varepsilon}^{0}$ is the midsurface strain and
$\boldsymbol{\kappa}$ is the curvature vector. Transverse shear strains are
included in the six-component material strain state. The layer stresses are
transformed to the shell local frame before the element force and stiffness
are assembled.

The element uses the layer thicknesses to define the total thickness:

```{math}
t = \sum_{i=1}^{n} t_i.
```

For a single-layer definition without `layers`, the implementation uses a
default thickness of `1.0` and a default orientation of $0^\circ$.

## Input parameters

The element block must contain `type = "SLS"` and a `material` block. A
single material can be used for all layers, or a `MultiMaterial` block can be
used when different layers refer to different material definitions.

### Single material

This is the compact form used by
[sls_cantilever01.pro](../../examples/elements/sls/sls_cantilever01.pro):

```text
SLSElem =
{
  type = "SLS";

  material =
  {
    type = "Isotropic";
    E    = 1.0e6;
    nu   = 0.0;
    rho  = 1.11e3;
  };
};
```

| Parameter | Description | Optional | Default value |
| --- | --- | --- | --- |
| `type` | Must be `"SLS"`. | No | — |
| `material` | Material definition used by the element. | No | — |
| `theta` | Material orientation in degrees. | Yes | $0^\circ$ |
| `E` | Young's modulus for an isotropic material. | Depends on material model | — |
| `nu` | Poisson's ratio for an isotropic material. | Depends on material model | — |
| `rho` | Density, used for the mass matrix. | Depends on material model | — |
| `incremental` | Selects the incremental material formulation when supported. | Yes | `false` |

The SLS implementation does not use a top-level `thickness` property in this
form. If no `layers` list is supplied, it creates one layer with thickness
`1.0`. To specify a physical thickness, use the layered form below.

### Layered material

Define the layer names in `layers`. Each layer specifies its thickness,
orientation, and material identifier. This pattern is used by
[sls_cantilever03.pro](../../examples/elements/sls/sls_cantilever03.pro).

```text
SLSElem =
{
  type   = "SLS";
  layers = ["c0", "c90", "c0"];

  c0  = { thickness = 1.0; theta = 0.0;  material = 0; };
  c90 = { thickness = 1.0; theta = 90.0; material = 0; };

  material =
  {
    type = "TransverseIsotropic";
    E1   = 1.0e6;
    E2   = 1.0e5;
    nu12 = 0.25;
    G12  = 1.0e5;
    rho  = 1.23e4;
  };
};
```

The element-level parameters are:

| Parameter | Description | Optional | Default value |
| --- | --- | --- | --- |
| `type` | Must be `"SLS"`. | No | — |
| `layers` | Ordered list of layer-block names through the thickness. | No | — |
| `material` | One material block, or a `MultiMaterial` block. | No | — |

Each layer block contains:

| Parameter | Description | Optional | Default value |
| --- | --- | --- | --- |
| `thickness` | Layer thickness. | No | — |
| `theta` | Layer orientation in degrees. | No in layered form | — |
| `material` | Material index or material name. | No | — |

For a single `material` block, the material properties depend on its type:

| Parameter | Description | Optional | Default value |
| --- | --- | --- | --- |
| `type` | Material model, for example `"Isotropic"` or `"TransverseIsotropic"`. | No | — |
| `E` | Isotropic Young's modulus. | Required for `Isotropic` | — |
| `nu` | Isotropic Poisson's ratio. | Required for `Isotropic` | — |
| `E1`, `E2` | Principal Young's moduli. | Required for `TransverseIsotropic` | — |
| `nu12` | Principal Poisson's ratio. | Required for `TransverseIsotropic` | — |
| `G12` | In-plane shear modulus. | Yes | $E_1/[2(1+\nu_{12})]$ |
| `G13` | 1–3 shear modulus. | Yes | `G12` |
| `G23` | 2–3 shear modulus. | Yes | `G12` |
| `rho` | Material density. | No for dynamic analysis | — |
| `incremental` | Uses the incremental formulation when supported. | Yes | `false` |

### Multiple materials

Use `type = "MultiMaterial"` when different layers use different material
blocks. The layer's `material` value must match one of the names in the
`materials` list. This is demonstrated by
[sls_cantilever04.pro](../../examples/elements/sls/sls_cantilever04.pro).

```text
material =
{
  type      = "MultiMaterial";
  materials = ["AA", "BB"];

  AA = { type = "TransverseIsotropic"; E1 = 1.0e6; E2 = 1.0e5;
         nu12 = 0.25; G12 = 1.0e5; rho = 1.23e4; };
  BB = { type = "TransverseIsotropic"; E1 = 1.0e6; E2 = 1.0e5;
         nu12 = 0.25; G12 = 1.0e5; rho = 2.23e4; };
};
```

## Output and examples

Stress output is generated by the selected material model. For the standard
six-component solid stress state, the labels are typically:

```text
S11, S22, S33, S23, S13, S12
```

The examples in `examples/elements/sls` cover the main use cases:

- [Single-material nonlinear cantilever](../../examples/elements/sls/sls_cantilever01.pro)
  demonstrates the basic SLS definition.
- [Oriented single-layer materials](../../examples/elements/sls/sls_cantilever02.pro)
  demonstrates separate SLS element blocks with different orientations and
  densities.
- [Layered transverse-isotropic shell](../../examples/elements/sls/sls_cantilever03.pro)
  demonstrates a three-layer definition using one material block.
- [Multi-material laminate](../../examples/elements/sls/sls_cantilever04.pro)
  demonstrates layers that refer to different material blocks.
- [Dynamic SLS cantilever](../../examples/elements/sls/sls_cantilever_dyn.pro)
  demonstrates eigenvalue analysis and the use of density.

The examples select either a `LinearSolver`, `NonlinearSolver`, or
`DynEigSolver` according to the analysis type.
