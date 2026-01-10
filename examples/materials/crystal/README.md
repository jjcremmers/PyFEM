# Crystal Plasticity Examples

This directory contains example input files for testing the `Crystal` material model in PyFEM.

## Model Overview

The `Crystal` material model implements rate-dependent crystal plasticity based on:
- Schmid's law for crystallographic slip
- Power-law viscoplasticity
- Self- and latent-hardening
- Multiple slip system sets
- Crystal orientation effects

## Examples

### 1. simple_tension.pro

**Description:** Basic uniaxial tension test of a single crystal aluminum specimen

**Features:**
- Single slip system set {110}<111>
- Identity crystal orientation
- Simple 2×2 mesh of quad elements
- Force-controlled loading

**Key Parameters:**
- Material: Aluminum-like properties
- Elastic: c11=108 GPa, c12=62 GPa, c44=28 GPa
- Plasticity: tau0=16 MPa, taus=148 MPa
- Rate sensitivity: m=20

**Usage:**
```bash
pyfem simple_tension.pro
```

**Expected Behavior:**
- Elastic response followed by plastic deformation
- Slip system activation depending on Schmid factors
- Rate-dependent stress-strain curve

---

### 2. oriented_shear.pro

**Description:** Shear test with 45-degree rotated crystal orientation

**Features:**
- Multiple slip system sets (3 sets covering {110}<111>)
- 45-degree crystal rotation about z-axis
- Larger mesh (6 elements)
- Shear loading configuration
- Tracks total accumulated plastic shear

**Key Parameters:**
- Material: Copper-like properties
- Elastic: c11=168 GPa, c12=121 GPa, c44=75 GPa
- Plasticity: tau0=60 MPa, taus=180 MPa
- Rate sensitivity: m=30
- Latent hardening: qab=1.4

**Usage:**
```bash
pyfem oriented_shear.pro
```

**Expected Behavior:**
- Multiple slip systems activated
- Crystal orientation affects which systems slip first
- Interaction between slip systems through latent hardening
- Higher rate sensitivity leads to smoother response

---

## Material Parameters

### Elastic Properties

**Option 1: Cubic Crystal (recommended for crystals)**
- `c11`: Elastic constant (normal stiffness)
- `c12`: Elastic constant (coupling)
- `c44`: Elastic constant (shear stiffness)

**Option 2: Isotropic (fallback)**
- `E`: Young's modulus
- `nu`: Poisson's ratio

### Slip Systems

- `nsets`: Number of slip system sets (1-3)
- `plane1`, `plane2`, `plane3`: Slip plane normals (Miller indices)
- `dir1`, `dir2`, `dir3`: Slip directions (Miller indices)

Common slip systems:
- FCC metals: {111}<110> or {110}<111>
- BCC metals: {110}<111>, {112}<111>, {123}<111>

### Crystal Orientation

- `orient1_local`, `orient1_global`: First orientation vector (in local and global coords)
- `orient2_local`, `orient2_global`: Second orientation vector

For identity orientation, use:
```
orient1_local  = [1.0, 0.0, 0.0];
orient1_global = [1.0, 0.0, 0.0];
orient2_local  = [0.0, 1.0, 0.0];
orient2_global = [0.0, 1.0, 0.0];
```

### Viscoplasticity

- `gamma0`: Reference shear strain rate (typically 0.001 s⁻¹)
- `m`: Rate sensitivity exponent (higher = more rate-dependent)
  - m=5-10: low rate sensitivity
  - m=20-30: moderate rate sensitivity
  - m>50: approaches rate-independent limit

### Hardening

- `tau0`: Initial critical resolved shear stress (CRSS)
- `h0`: Initial hardening modulus
- `taus`: Saturation stress
- `a`: Hardening exponent (controls hardening rate decay)
- `qab`: Latent hardening ratio (qab=1.0 → self-hardening only, qab>1.0 → latent hardening)

**Hardening Law (Bassani-Wu):**
```
h_αβ = h0 * q_αβ * (1 - tau_c/tau_s)^a
```

### Integration Parameters

- `theta`: Implicit integration parameter (0.5 = mid-point, recommended)
- `nlgeom`: Use finite deformation (false for small strain)
- `maxiter`: Maximum Newton-Raphson iterations (20-25 typical)
- `tol`: Convergence tolerance (1e-6 typical)

---

## Solver Settings

**Important for rate-dependent materials:**
- Use `NonlinearSolver` (not `LinearSolver`)
- Set appropriate `dtime` (time increment)
  - Smaller dtime → more accurate rate effects
  - Typical: 0.01 to 0.1 depending on loading rate
- May need smaller convergence tolerance for complex slip patterns

---

## Output

The model provides the following output quantities:

**Stress Components:**
- S11, S22, S33: Normal stresses
- S23, S13, S12: Shear stresses

**Plastic Deformation:**
- GammaTotal: Total accumulated plastic shear (sum over all systems)
- Gamma1, Gamma2, ...: Shear strain in individual slip systems
- Tau1, Tau2, ...: Resolved shear stress in individual slip systems

---

## Tips for Creating Custom Examples

1. **Start simple:** Begin with 1 slip system set and identity orientation
2. **Check Schmid factors:** Verify slip systems are oriented to activate under your loading
3. **Tune dtime:** If convergence fails, try smaller time increments
4. **Monitor slip activity:** Use GraphWriter to track Gamma1, Tau1, etc.
5. **Validate elasticity:** Check elastic behavior matches c11, c12, c44 before adding plasticity

---

## Material-Specific Suggestions

### Aluminum (FCC)
```
c11 = 108e3;  c12 = 62e3;  c44 = 28e3;
tau0 = 16;    taus = 148;   h0 = 180;
```

### Copper (FCC)
```
c11 = 168e3;  c12 = 121e3;  c44 = 75e3;
tau0 = 60;    taus = 180;   h0 = 200;
```

### Nickel (FCC)
```
c11 = 247e3;  c12 = 147e3;  c44 = 125e3;
tau0 = 80;    taus = 200;   h0 = 300;
```

### Iron/Steel (BCC)
```
c11 = 230e3;  c12 = 135e3;  c44 = 117e3;
tau0 = 100;   taus = 300;   h0 = 400;
```

---

## References

- Peirce, D., Asaro, R. J., & Needleman, A. (1983). "Material rate dependence and localized deformation in crystalline solids." Acta Metallurgica, 31(12), 1951-1976.
- Bassani, J. L., & Wu, T. Y. (1991). "Latent hardening in single crystals. II. Analytical characterization and predictions." Proceedings of the Royal Society of London A, 435(1893), 21-41.
- Kysar, J. W. (1997). "Addendum to 'A user-material subroutine incorporating single crystal plasticity in the ABAQUS finite element program'."

---

## Troubleshooting

**Problem:** Convergence failure
- Solution: Reduce dtime, increase maxiter, or adjust convergence tolerance

**Problem:** No plastic deformation
- Solution: Check Schmid factors - slip systems may not be favorably oriented

**Problem:** Excessive plastic deformation
- Solution: Increase tau0 or check loading magnitude

**Problem:** Negative determinant
- Solution: Use nlgeom=false for small strain, or check element distortion
