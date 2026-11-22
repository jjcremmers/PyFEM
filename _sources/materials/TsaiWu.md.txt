# TsaiWu

`TsaiWu` is a quadratic failure criterion for orthotropic and fiber-reinforced
composite materials. It accounts for interaction between the normal stresses
and distinguishes tensile and compressive strengths in the principal material
directions.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Failure index

For a plane-stress vector $[\sigma_1,\sigma_2,\tau_{12}]$, the implementation
evaluates

```{math}
FI = F_1\sigma_1 + F_2\sigma_2 + F_{11}\sigma_1^2
 + F_{22}\sigma_2^2 + F_{66}\tau_{12}^2
 + 2F_{12}\sigma_1\sigma_2.
```

The returned strength ratio is $\sqrt{FI}$ for non-negative `FI`; a negative
value is returned when the computed index is negative. A value of 1 indicates
the failure surface.

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"TsaiWu"`. | String | - |
| `Xt`, `Xc` | Tensile and compressive strengths in direction 1. | Float | - |
| `Yt`, `Yc` | Tensile and compressive strengths in direction 2. | Float | - |
| `S` | In-plane shear strength. | Float | - |
| `F12` | Optional normal-stress interaction coefficient. | Float | Optional |

If `F12` is omitted, the implementation uses
`-0.5*sqrt(F11*F22)`. Stresses must be supplied in the material coordinate
system.
