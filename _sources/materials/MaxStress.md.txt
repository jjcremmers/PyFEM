# MaxStress

`MaxStress` is a simple composite failure criterion. It compares the normal
and shear stresses in the material directions with separate tensile,
compressive, and shear strengths.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Failure index

The criterion evaluates directions 1 and 2 separately, using `Xt`/`Xc` and
`Yt`/`Yc` for tension and compression. The absolute in-plane shear stress is
compared with `S`. The largest normalized value is returned as the failure
index; values greater than or equal to 1 indicate failure.

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"MaxStress"`. | String | Required |
| `Xt`, `Xc` | Tensile and compressive strengths in direction 1. | Float | Required |
| `Yt`, `Yc` | Tensile and compressive strengths in direction 2. | Float | Required |
| `S` | Shear strength. | Float | Required |

The stress components must be expressed in the material coordinate system.
