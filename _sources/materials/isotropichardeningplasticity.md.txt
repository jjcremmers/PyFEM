# IsotropicHardeningPlasticity

The `IsotropicHardeningPlasticity` model describes J2 plasticity with isotropic hardening. Requires elastic properties and a yield stress, with a hardening law defined via the Hardening utility.

## Overview
- **Material type:** `IsotropicHardeningPlasticity`

## Parameters
### Mandatory
- `type`: Must be set to `"IsotropicHardeningPlasticity"`
- `E`: Young's modulus
- `nu`: Poisson's ratio
- `syield`: Initial yield stress
