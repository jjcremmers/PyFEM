# Dummy

`Dummy` is a simple linear traction law intended primarily for testing
interface elements. It applies the same scalar stiffness `D` to the normal and
tangential components of the interface displacement jump. It is useful for
verification and basic interface experiments, but it does not represent a
physical cohesive failure law.

```{warning}
This section was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used. It is documentation support, not a substitute for verification of the
underlying equations or code.
```

## Overview

- **Material type:** `Dummy`
- **Application:** interface-element testing
- **Supported dimensions:** 2D and 3D interface elements
- **Traction output:** `Tn`, `Ts` in 2D; `Tn`, `Ts1`, `Ts2` in 3D

The material uses the interface displacement jump as its deformation vector
and computes traction as

```{math}
\boldsymbol{t} = \boldsymbol{H}\,\boldsymbol{\delta},
\qquad
\boldsymbol{H} = D\boldsymbol{I}.
```

The tangent is constant and diagonal. No damage, softening, or failure
evolution is included.

## Input parameters

```text
material =
{
  type = "Dummy";
  D    = 1.0e5;
};
```

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"Dummy"`. | String | Required |
| `D` | Interface stiffness applied to each traction component. | Float | Required |

The material must be assigned to an [`Interface` element](../elements/interface.md).

## Examples

- [Traction oscillations example](../../examples/ch13/TractionOscillations.pro)
- [Interface element documentation](../elements/interface.md)
