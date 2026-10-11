# PowerLawModeI

`PowerLawModeI` is a mode-I cohesive material for interface elements. It
relates the normal opening displacement to the normal traction using a smooth
exponential law. The response rises from zero to a peak traction and then
softens as the interface separates. The fracture energy `Gc` controls the area
under the traction–separation curve, while `Tult` controls its peak strength.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used. It is documentation support, not a substitute for verification of the
underlying equations or code.
```

## Overview

- **Material type:** `PowerLawModeI`
- **Application:** mode-I opening and debonding
- **Traction output:** `Tn`, `Ts`
- **Active response:** normal component only; tangential traction remains zero

The material must be assigned to an [`Interface` element](../elements/interface.md).
The normal opening is read from the first component of the interface
deformation vector. The characteristic opening is

```{math}
\delta_n = \frac{G_c}{e\,T_{ult}},
```

and the implemented normal traction law is

```{math}
T_n(\delta) = \frac{G_c}{\delta_n^2}
\exp\left(-\frac{\delta}{\delta_n}\right)\delta.
```

The returned tangent is the derivative of this expression with respect to the
normal opening. The implementation does not model tangential separation.

## Input parameters

```text
material =
{
  type = "PowerLawModeI";
  Gc   = 1.0;
  Tult = 10.0;
};
```

| Parameter | Description | Type |
| --- | --- | --- |
| `type` | Must be `"PowerLawModeI"`. | String |
| `Gc` | Fracture energy. | Float |
| `Tult` | Ultimate or peak normal traction. | Float |

## Examples

- [Peel test](../../examples/ch13/PeelTest.pro)
- [Peel test with 80 elements](../../examples/ch13/PeelTest80.pro)
- [Peel test with 100 elements](../../examples/ch13/PeelTest100.pro)
- [Interface element documentation](../elements/interface.md)
