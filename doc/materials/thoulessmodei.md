# ThoulessModeI

`ThoulessModeI` is a piecewise-linear mode-I cohesive material for interface
elements. It describes an initial elastic response, a constant-traction region
at the peak strength, and a linear softening branch leading to complete
decohesion. The characteristic openings are derived from the fracture energy,
peak traction, and the ratios `d1d3` and `d2d3`.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used. It is documentation support, not a substitute for verification of the
underlying equations or code.
```

## Overview

- **Material type:** `ThoulessModeI`
- **Application:** mode-I opening and debonding
- **Traction output:** `Tn`, `Ts`
- **Active response:** normal component only; tangential traction remains zero

The material must be assigned to an [`Interface` element](../elements/interface.md).
The normal opening is read from the first component of the interface
deformation vector. The derived openings are

```{math}
\delta_3 = \frac{2G_c}{(-r_1+r_2+1)T_{ult}},
\qquad
\delta_1 = r_1\delta_3,
\qquad
\delta_2 = r_2\delta_3,
```

where $r_1$ is `d1d3` and $r_2$ is `d2d3`. The response is linear up to
$\delta_1$, constant at `Tult` between $\delta_1$ and $\delta_2$, and linearly
softens to zero at $\delta_3$. Openings beyond $\delta_3$ produce zero normal
traction.

## Input parameters

```text
material =
{
  type = "ThoulessModeI";
  Gc   = 1.0;
  Tult = 10.0;
  d1d3 = 0.1;
  d2d3 = 0.5;
};
```

| Parameter | Description | Type |
| --- | --- | --- |
| `type` | Must be `"ThoulessModeI"`. | String |
| `Gc` | Fracture energy. | Float |
| `Tult` | Peak normal traction. | Float |
| `d1d3` | Ratio \(\delta_1/\delta_3\). | Float |
| `d2d3` | Ratio \(\delta_2/\delta_3\). | Float |

## Examples

There is currently no dedicated `ThoulessModeI` example in the repository. It
can be configured in the same way as the other mode-I materials in the
[peel-test examples](../../examples/ch13/PeelTest.pro), after verifying that
the selected parameters and interface convention are appropriate.
