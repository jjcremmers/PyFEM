# Interface

The `Interface` element models cohesive or contact interfaces between continuum parts. It is used to simulate delamination, peeling, and traction-separation behaviors using appropriate interface constitutive laws.

## Overview
- **Element type:** `Interface`
- **Supports:**
  - Traction–separation laws via material models (e.g., `XuNeedleman`)
  - 2D and 3D interface formulations
  - Nonlinear analysis for softening and debonding

## Parameters
### Mandatory
- `type`: Must be set to `"Interface"`
- `material`: Interface material model block defining cohesive law and parameters.

### Optional
(See documentation for details)
