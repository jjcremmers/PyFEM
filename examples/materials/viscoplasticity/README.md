# Viscoplastic Material Example

## Overview

This example demonstrates rate-dependent plastic behavior using the `ViscoPlasticity` material model based on the Perzyna overstress formulation.

## Problem Description

A notched bar in tension is simulated where:
- The bar has a geometric notch creating stress concentration
- The material exhibits rate-dependent yielding
- Higher loading rates result in higher apparent yield stress
- Plastic deformation localizes at the notch

## Material Model

The `ViscoPlasticity` model uses Perzyna overstress theory:

```
d(εᵖ)/dt = γ × <(σ_vm - σ_y)/σ_y>ⁿ × N
```

where:
- γ: fluidity parameter (controls rate sensitivity)
- σ_vm: von Mises stress
- σ_y: current yield stress (with hardening)
- n: rate sensitivity exponent
- N: flow direction (normal to yield surface)

### Material Parameters

- **E = 200000.0**: Young's modulus (MPa) - typical for steel
- **nu = 0.3**: Poisson's ratio
- **syield = 250.0**: Initial yield stress (MPa)
- **hard = 1000.0**: Linear hardening modulus (MPa)
- **gamma = 0.001**: Fluidity parameter (1/MPa/s)
- **n = 1.0**: Rate sensitivity exponent (linear viscosity)

## Expected Behavior

1. **Elastic loading**: Initially, the material responds elastically
2. **Rate-dependent yielding**: When stress exceeds yield, plastic flow rate depends on:
   - Magnitude of overstress (σ - σ_y)
   - Loading rate (through time step and gamma)
3. **Stress concentration**: Maximum stress and plastic strain occur at the notch
4. **Hardening**: Yield stress increases with plastic deformation
5. **Smooth transition**: No sharp yield point - gradual onset of plasticity

## Physical Interpretation

### Rate Dependency
- **Fast loading** (small dtime): Material appears stiffer, higher peak stress
- **Slow loading** (large dtime): More plastic flow, lower peak stress
- **gamma → 0**: Approaches rate-independent plasticity
- **gamma → ∞**: Approaches viscous flow

### Applications
This model is relevant for:
- **Metals at high strain rates**: Impact, crash simulations
- **Polymers**: Highly rate-dependent behavior
- **High temperature metals**: Creep and viscoplasticity
- **Numerical regularization**: Avoiding mesh-dependent localization

## Running the Example

```bash
pyfem bar_tension.pro
```

## Output

The GraphWriter module tracks:
- **time**: Simulation time
- **disp**: Horizontal displacement at the right end (node 17)
- **load**: Total reaction force at loaded nodes
- **plastic**: Equivalent plastic strain at the notch (node 8)

## Geometry Details

The notched bar has:
- Length: 8.0 units
- Width: 0.5 units (0.45 at notch)
- Thickness: 0.2 units
- Notch location: Center (x = 3.5-4.0)
- Notch depth: 0.05 units (10% width reduction)

## Variations to Try

1. **Different loading rates**:
   ```
   dtime = 0.001;  # Fast loading - higher peak stress
   dtime = 0.1;    # Slow loading - more creep
   ```

2. **Rate sensitivity exponent**:
   ```
   n = 2.0;  # More rate-sensitive at high overstress
   n = 0.5;  # Less rate-sensitive
   ```

3. **Fluidity parameter**:
   ```
   gamma = 0.0001;  # Less viscous, closer to rate-independent
   gamma = 0.01;    # More viscous, significant rate effects
   ```

4. **Compare with rate-independent**:
   - Set gamma very small (e.g., 1e-10)
   - Should match `IsotropicKinematicHardening` results

## Understanding the Results

**Load-displacement curve**:
- Initially linear (elastic)
- Smooth transition to plastic flow (no sharp yield point)
- Strain hardening causes rising load
- Rate effects visible when comparing different dtime values

**Plastic strain distribution**:
- Localizes at the notch (stress concentration)
- Gradual spread with continued loading
- Regularized by viscosity (gamma parameter)
