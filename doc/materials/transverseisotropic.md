# TransverseIsotropic

`TransverseIsotropic` is a three-dimensional linear-elastic material model for
unidirectional composites with one preferred material direction. The material
uses `E1` in the preferred direction, `E2` in the transverse directions, and
the shear modulus `G12` for the shear terms.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used. It is documentation support, not a substitute for verification of the
underlying equations or code.
```

## Overview

- **Material type:** `TransverseIsotropic`
- **Kinematics:** three-dimensional linear elasticity
- **Stress output:** `S11`, `S22`, `S33`, `S23`, `S13`, `S12`
- **Formulation:** total strain by default

The stress and strain vectors use the Voigt ordering

```text
[S11, S22, S33, S23, S13, S12]
```

The preferred material direction is component 1. Components 2 and 3 define
the transverse plane. The implementation returns a constant six-by-six
constitutive matrix as the material tangent.

## Constitutive response

The stress is computed from the strain vector as

```{math}
\boldsymbol{\sigma} = \boldsymbol{H}\,\boldsymbol{\varepsilon}.
```

The normal-stress block is determined by `E1`, `E2`, and `nu12`. The three
shear entries of the implemented matrix are set to `G12`:

```{math}
H_{44} = H_{55} = H_{66} = G_{12}.
```

The material assumes that the material axes are aligned with the element
coordinate system. Verify this assumption before using the model for rotated
laminae or general composite layups.

## Input parameters

```text
material =
{
  type = "TransverseIsotropic";
  E1   = 135.0;
  E2   = 10.0;
  nu12 = 0.3;
  G12  = 5.0;
};
```

| Parameter | Description | Type |
| --- | --- | --- |
| `type` | Must be `"TransverseIsotropic"`. | String |
| `E1` | Young's modulus in the preferred material direction. | Float |
| `E2` | Young's modulus in the transverse directions. | Float |
| `nu12` | Major Poisson's ratio. | Float |
| `G12` | Shear modulus used by the constitutive matrix. | Float |

Use a consistent unit system for all elastic parameters and the resulting
stresses.

## Examples

- [SLS cantilever with transverse-isotropic materials](../../examples/elements/sls/sls_cantilever02.pro)
- [SLS cantilever with multiple material regions](../../examples/elements/sls/sls_cantilever04.pro)
- [Composite cantilever exercise](../../exercises/composites/cantilever_sls_carb.pro)
- [Multiple-material composite exercise](../../exercises/composites/cantilever_sls_multiple.pro)
