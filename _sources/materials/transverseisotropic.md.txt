# TransverseIsotropic

The `TransverseIsotropic` material model describes linear elastic behavior of unidirectional composites with one preferred fiber direction. Outputs 3D stress components: S11, S22, S33, S23, S13, S12.

## Overview
- **Material type:** `TransverseIsotropic`

## Parameters
### Mandatory
- `type`: Must be set to `"TransverseIsotropic"`
- `E1`: Young's modulus along fiber direction
- `E2`: Young's modulus transverse to fiber
- `nu12`: Major Poisson's ratio
- `G12`: In-plane shear modulus
