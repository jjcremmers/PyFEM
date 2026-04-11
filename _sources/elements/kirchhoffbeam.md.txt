# KirchhoffBeam

The `KirchhoffBeam` element models slender beams using the Kirchhoff (Euler) beam theory, suitable for frames and slender structures where shear deformation can be neglected.

## Overview
- **Element type:** `KirchhoffBeam`
- **Supports:**
  - Frame analysis with bending and axial deformation
  - Suitable for slender beams (no shear deformation)
  - 2D frame applications

## Parameters
### Mandatory
- `type`: Must be set to `"KirchhoffBeam"`
- `material`: Material block with elastic properties (e.g., `E`, `rho`) and geometric properties via mesh/section data.

### Optional
(See documentation for details)
