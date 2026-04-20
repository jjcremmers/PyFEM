# Reissner-Mindlin Shell

The `ReissnerMindlinShell` element implements a 4-node shell formulation with transverse shear deformation, 6 degrees of freedom per node, and geometrically nonlinear director-based kinematics.

## Overview
- **Element type:** `ReissnerMindlinShell`
- **Degrees of freedom per node:** `u`, `v`, `w`, `rx`, `ry`, `rz`
- **Current scope:**
  - 4-node quadrilateral midsurfaces only
  - Initially curved or flat midsurfaces
  - Layer-wise shell integration through the thickness using the same laminate setup as `Plate`
  - Geometrically nonlinear midsurface/director kinematics
  - A small drilling stabilization for the `rz` degree of freedom

## Parameters
### Mandatory
- `type`: Must be set to `"ReissnerMindlinShell"`
- `material`: Single-material block or laminate setup as used by `Plate`
- `thickness`: Shell thickness for single-layer setups

### Optional
- `layers`: Laminate layer list
- `materials`: Laminate material list
- `shearCorrection`: Shear correction factor
- `tangentPerturbation`: Finite-difference perturbation used for the numerical tangent
- `drillingScale`: Scaling factor for the drilling stabilization

## Notes
- The current implementation supports curved Quad4 midsurfaces, but it is still restricted to that interpolation.
- The tangent stiffness is assembled from material and geometric contributions using the current shell kinematics, with the director Jacobian evaluated locally.
