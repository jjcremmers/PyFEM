# SandwichCore material

`SandwichCore` is a linear elastic material model for sandwich-panel cores.
It represents a core that is comparatively soft in the in-plane directions
and stiffer in the through-thickness direction. The model is useful for
idealized foam, honeycomb, and other lightweight core materials when a simple
orthotropic-like constitutive response is sufficient.

The material is defined by a through-thickness Young's modulus `E3`, a base
shear modulus `G`, and an in-plane stiffness scaling factor. The default factor
reduces the in-plane normal stiffness to `0.1%` of `E3`.

`SandwichCore` is a linear elastic material model. It does not include
plasticity, damage, failure, or history variables. Use a material with the
required nonlinear behavior if the core response must evolve during the
analysis.


## Implementation

The implementation is in `pyfem/materials/SandwichCore.py`. The material uses
the six-component Voigt ordering

```text
[S11, S22, S33, S23, S13, S12]
```

and computes the stress from the engineering strain vector as

```{math}
\boldsymbol{\sigma} = \boldsymbol{H}\,\boldsymbol{\varepsilon}.
```

The constitutive matrix is diagonal:

```{math}
\boldsymbol{H} =
\begin{bmatrix}
fE_3 & 0 & 0 & 0 & 0 & 0 \\
0 & fE_3 & 0 & 0 & 0 & 0 \\
0 & 0 & E_3 & 0 & 0 & 0 \\
0 & 0 & 0 & f(G_{13}+G_{23})/2 & 0 & 0 \\
0 & 0 & 0 & 0 & G_{23} & 0 \\
0 & 0 & 0 & 0 & 0 & G_{13}
\end{bmatrix},
```

where $f$ is the `factor` parameter. The default is $f=0.001$. The material
returns this matrix as its consistent tangent and reports the computed stress
components as output data.

:::{warning}
The shear terms follow the component ordering and assignments implemented in
`SandwichCore.py`. Verify the resulting material response against the chosen
element's strain and stress convention before using non-default `G13` and
`G23` values.
:::

## Input parameters

The material is defined inside an element material block. A minimal
configuration is:

```text
material =
{
  type = "SandwichCore";
  E3   = 100.0;
  G    = 10.0;
};
```

The available parameters are:

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"SandwichCore"`. | String | — |
| `E3` | Through-thickness Young's modulus. | Float | — |
| `G` | Base shear modulus. | Float | — |
| `factor` | Scaling factor for the in-plane normal stiffness and one shear contribution. | Float | Optional, default = 0.001 |
| `G13` | Shear modulus supplied to the constitutive matrix. | Float | Optional, default = `G` |
| `G23` | Shear modulus supplied to the constitutive matrix. | Float | Optional, default = `G` |

All material parameters must use a consistent unit system. For example, if
`E3` and `G` are specified in force per area, the resulting stresses have the
same units.

## Examples

There is currently no dedicated `SandwichCore` example in the repository. The
material can be added to a supported continuum element as follows:

```text
CoreElem =
{
  type = "FiniteStrainContinuum";

  material =
  {
    type   = "SandwichCore";
    E3     = 100.0;
    G      = 10.0;
    factor = 0.001;
  };
};
```

For a more detailed example, provide a mesh and boundary conditions for the
selected continuum element, then choose an appropriate linear or nonlinear
solver. The material itself is linear, but it can be used within a nonlinear
analysis when the surrounding model or geometry requires it.
