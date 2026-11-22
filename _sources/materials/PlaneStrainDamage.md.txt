# PlaneStrainDamage

`PlaneStrainDamage` is a two-dimensional plane-strain damage model based on
linear elasticity and an equivalent-strain history variable. Damage grows
irreversibly once the equivalent strain exceeds the stored maximum and reduces
the effective stress and tangent stiffness.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Overview

The model returns `S11`, `S22`, `S12`, and `damage`. The damage variable ranges
from 0 for the undamaged state to 1 for complete degradation. The response is
intended for plane-strain continuum elements.

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"PlaneStrainDamage"`. | String | Required |
| `E` | Young's modulus. | Float | Required |
| `nu` | Poisson's ratio. | Float | Required |
| `k` | Damage-model parameter controlling the equivalent-strain measure. | Float | Required |
| `kappa0` | Damage-initiation equivalent strain. | Float | Required |
| `kappac` | Equivalent strain at complete damage. | Float | Required |

The model stores the history variable `kappa`; unloading does not reduce the
maximum equivalent strain or the accumulated damage.

## Examples

- [Continuum damage example](../../examples/ch06/ContDamExample.pro)
- [Micromodel damage example](../../examples/materials/micromodel/micro_rve.pro)
