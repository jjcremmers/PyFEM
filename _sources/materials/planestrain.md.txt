# PlaneStrain

The `PlaneStrain` material model implements linear elastic behavior in 2D plane strain. Use for thick structures or long domains where out-of-plane strain is negligible. Outputs: S11, S22, S12.

## Overview
- **Material type:** `PlaneStrain`

## Parameters
### Mandatory
- `type`: Must be set to `"PlaneStrain"`
- `E`: Young's modulus
- `nu`: Poisson's ratio
