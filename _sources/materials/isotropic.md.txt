# Isotropic

The `Isotropic` material model implements linear elastic, isotropic behavior for 3D solids. It can be used with continuum elements in both small strain and finite strain analyses where the constitutive response is governed by Young's modulus and Poisson's ratio. Density can be provided for dynamic problems. TEST

## Overview
- **Material type:** `Isotropic`
- Supports 3D linear elastic isotropic behavior
- Compatible with small strain and finite strain kinematics (via element choice)
- Stress output: S11, S22, S33, S23, S13, S12

## Parameters
### Mandatory
- `type`: Must be set to `"Isotropic"`
- `E`: Young's modulus
- `nu`: Poisson's ratio
