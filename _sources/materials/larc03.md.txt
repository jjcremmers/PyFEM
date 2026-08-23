# Larc03

`Larc03` implements the LARC03 failure criterion for unidirectional composite
plies under quasi-2D stress states. It distinguishes between matrix tension,
matrix compression, fiber tension, and fiber compression, and uses a
misaligned-frame check for fiber compression.

## Overview
- **Failure type:** `Larc03`
- **Configured through:** `failureType` in the material block
- **Typical use:** ply-level failure checks for orthotropic composite materials

## Parameters
### Mandatory elastic constants
- `E1`: Young's modulus along the fiber direction
- `E2`: Young's modulus transverse to the fiber direction
- `G12`: In-plane shear modulus
- `nu12`: Major Poisson ratio

### Mandatory strengths
- `Xt` or `XT` or `F1t`: fiber tensile strength
- `Xc` or `XC` or `F1c`: fiber compressive strength
- `Yt` or `YT` or `F2t`: transverse tensile strength
- `Yc` or `YC` or `F2c`: transverse compressive strength
- `S` or `SL` or `Fs`: in-plane shear strength

### In-situ data
Provide one of the following two sets:

- Direct in-situ strengths:
  `YT_is`, `SL_is`
- Fracture-data-based input:
  `GIc_L`, `GIIc_L`, `ply_thickness`

### Optional parameters
- `etaL`: longitudinal friction parameter used in compression
- `g`: mixed-mode weighting factor for matrix tension
- `alpha0`: fracture-plane angle in radians
- `alpha0_deg`: fracture-plane angle in degrees when `alpha0` is not given

## Example

```text
material = {
  type        = "TransverseIsotropic";
  failureType = "Larc03";

  E1   = 135000.0;
  E2   = 10000.0;
  G12  = 5000.0;
  nu12 = 0.30;

  Xt = 1500.0;
  Xc = 1200.0;
  Yt = 45.0;
  Yc = 200.0;
  S  = 90.0;

  GIc_L        = 0.30;
  GIIc_L       = 0.90;
  ply_thickness = 0.125;
  alpha0_deg   = 53.0;
};
```

## Notes
- The failure criterion expects stresses in the material coordinate system.
- For 2D states it uses `[sigma11, sigma22, tau12]`.
- For 3D states it reads `tau12` from the sixth stress component.
- The returned value is a failure index; values greater than or equal to `1.0`
  indicate failure initiation.