# PlaneStress

The `PlaneStress` material model implements linear elastic behavior in 2D plane stress. Use for thin structures where out-of-plane stress is negligible. Outputs stress components: S11, S22, S12.

## Overview
- **Material type:** `PlaneStress`

## Parameters
### Mandatory
- `type`: Must be set to `"PlaneStress"`
- `E`: Young's modulus
- `nu`: Poisson's ratio
