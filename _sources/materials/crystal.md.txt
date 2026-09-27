# Crystal

The `Crystal` material model implements a small-strain, rate-dependent crystal
plasticity formulation for single crystals. Plastic flow occurs through slip on
user-defined crystallographic systems, with optional crystal orientation,
power-law slip kinetics, and self/latent hardening.

## Overview
- **Material type:** `Crystal`
- **Typical use:** single-crystal or strongly textured material response with slip-system-based plasticity
- **Kinematics:** small strain

## Parameters
### Elastic data
- `c11`, `c12`, `c44`: cubic-crystal elastic constants
- `E`, `nu`: isotropic elastic fallback when cubic constants are not provided

### Slip-system definition
- `nsets`: number of slip-system sets
- `plane1`, `plane2`, `plane3`: slip plane normals
- `dir1`, `dir2`, `dir3`: slip directions associated with each set

### Crystal orientation
- `orient1_local`, `orient1_global`: first reference vector in crystal and global coordinates
- `orient2_local`, `orient2_global`: second reference vector in crystal and global coordinates

### Slip kinetics and hardening
- `gamma0`: reference slip rate
- `m`: rate sensitivity exponent
- `tau0`: initial critical resolved shear stress
- `h0`: initial hardening modulus
- `taus`: saturation strength
- `a`: hardening exponent
- `qab`: latent-hardening ratio

### Integration controls
- `theta`: implicit integration parameter
- `maxiter`: local Newton iteration limit
- `tol`: local convergence tolerance

## Example

```text
material =
{
  type = "Crystal";

  c11 = 108.0e3;
  c12 = 62.0e3;
  c44 = 28.0e3;

  nsets  = 1;
  plane1 = [1.0, 1.0, 0.0];
  dir1   = [1.0, 1.0, 1.0];

  orient1_local  = [1.0, 0.0, 0.0];
  orient1_global = [1.0, 0.0, 0.0];
  orient2_local  = [0.0, 1.0, 0.0];
  orient2_global = [0.0, 1.0, 0.0];

  gamma0 = 0.001;
  m      = 20.0;
  tau0   = 16.0;
  h0     = 180.0;
  taus   = 148.0;
  a      = 2.25;
  qab    = 1.4;

  theta   = 0.5;
  maxiter = 20;
  tol     = 1.0e-6;
};
```

## Example Files
- `examples/materials/crystal/minimal_single_slip.pro`: minimal one-element smoke test
- `examples/materials/crystal/simple_tension.pro`: tensile loading example
- `examples/materials/crystal/oriented_shear.pro`: rotated-crystal shear example

## Notes
- The implementation works internally in 3D Voigt notation and maps back to 2D when used with 2D continuum elements.
- Output fields include stress components, accumulated slip, and slip-system shear variables for the first few systems.