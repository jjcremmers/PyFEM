# ThoulessModeI

`ThoulessModeI` defines a piecewise-linear mode-I traction–separation law with three characteristic displacement thresholds.

## Overview
- **Material type:** `ThoulessModeI`

## Parameters
### Mandatory
- `type`: Must be set to `"ThoulessModeI"`
- `Gc`: Fracture energy
- `Tult`: Peak traction
- `d1d3` / `d2d3`: Ratios defining the displacement thresholds

## Examples
- Cohesive interface modeling with piecewise traction law
