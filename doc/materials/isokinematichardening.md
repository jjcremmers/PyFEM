# IsotropicKinematicHardening

The `IsotropicKinematicHardening` model describes J2 plasticity with kinematic hardening (backstress evolution). Requires elastic properties, yield stress, and a kinematic hardening modulus.

## Overview
- **Material type:** `IsotropicKinematicHardening`

## Parameters
### Mandatory
- `type`: Must be set to `"IsotropicKinematicHardening"`
- `E`: Young's modulus
- `nu`: Poisson's ratio
- `syield`: Yield stress
- `hard`: Kinematic hardening modulus
