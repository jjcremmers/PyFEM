# TimoshenkoBeam

The `TimoshenkoBeam` element models beams including shear deformation effects, providing improved accuracy for moderately thick beams compared to Kirchhoff beam theory.

## Overview
- **Element type:** `TimoshenkoBeam`
- **Supports:**
  - Bending with shear deformation
  - Frame analysis where shear effects are non-negligible

## Parameters
### Mandatory
- `type`: Must be set to `"TimoshenkoBeam"`
- `material`: Material block with elastic properties (e.g., `E`, `G`, `rho`).

### Optional
(See documentation for details)
