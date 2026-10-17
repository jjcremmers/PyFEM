# IsotropicKinematicHardening

The `IsotropicKinematicHardening` model describes J2 plasticity with kinematic hardening (backstress evolution). Requires elastic properties, yield stress, and a kinematic hardening modulus.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Overview
- **Material type:** `IsotropicKinematicHardening`

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"IsotropicKinematicHardening"`. | String | Required |
| `E` | Young's modulus. | Float | Required |
| `nu` | Poisson's ratio. | Float | Required |
| `syield` | Initial yield stress. | Float | Required |
| `hard` | Kinematic hardening modulus. | Float | Required |
