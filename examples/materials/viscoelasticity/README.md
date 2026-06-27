# Viscoelastic Material Example

## Overview

This example demonstrates the time-dependent behavior of viscoelastic materials using the `ViscoElasticity` material model based on the generalized Maxwell model (Prony series representation).

## Problem Description

A simple single-element creep test is simulated where:
- A 3D cubic element is subjected to constant tensile load
- The displacement response is tracked over time
- The material exhibits time-dependent creep behavior

## Material Model

The `ViscoElasticity` model uses a generalized Maxwell representation:

```
E(t) = E_inf + Σ E_i × exp(-t/τ_i)
```

### Material Parameters

- **E = 1000.0**: Instantaneous Young's modulus (total stiffness at t=0)
- **nu = 0.3**: Poisson's ratio
- **Einf = 100.0**: Long-term equilibrium modulus (stiffness at t→∞)
- **nMaxwell = 3**: Number of Maxwell elements
- **relaxTimes = [0.1, 1.0, 10.0]**: Relaxation times for each Maxwell element
- **relaxModuli = [300.0, 300.0, 300.0]**: Relaxation moduli for each element

The sum: Einf + Σ E_i = 100 + 300 + 300 + 300 = 1000 = E ✓

## Expected Behavior

1. **Immediate elastic response**: Upon loading, the material responds with instantaneous modulus E
2. **Viscoelastic creep**: Over time, the material "relaxes" as the Maxwell elements deform
3. **Long-term behavior**: Eventually reaches equilibrium with effective modulus Einf
4. **Displacement increase**: The displacement gradually increases even under constant load

## Running the Example

```bash
pyfem creep_test.pro
```

## Output

The GraphWriter module tracks:
- **time**: Simulation time
- **disp**: Horizontal displacement at node 2
- **stress**: Stress component S11 at node 1

## Physical Interpretation

This example simulates creep behavior common in:
- Polymers at room temperature
- Metals at elevated temperatures
- Concrete under sustained loading
- Biological soft tissues

The multiple relaxation times capture the multi-scale nature of viscoelastic relaxation processes.

## Variations to Try

1. **Different loading rates**: Modify dtime to see rate-dependent stiffness
2. **More Maxwell elements**: Increase nMaxwell for smoother relaxation
3. **Different time scales**: Adjust relaxTimes to match specific material behavior
4. **Cyclic loading**: Add loading/unloading cycles to see hysteresis
