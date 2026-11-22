# IsotropicHardeningPlasticity

The `IsotropicHardeningPlasticity` model describes J2 plasticity with isotropic hardening. Requires elastic properties and a yield stress, with a hardening law defined via the Hardening utility.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Overview
- **Material type:** `IsotropicHardeningPlasticity`

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"IsotropicHardeningPlasticity"`. | String | Required |
| `E` | Young's modulus. | Float | Required |
| `nu` | Poisson's ratio. | Float | Required |
| `syield` | Initial yield stress. | Float | Required |
