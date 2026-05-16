# MultiMaterial

The `MultiMaterial` model groups multiple material models for layered or zone-based analyses. It delegates `getStress()` to the active sub-material specified by the element (e.g., via `iMat`).

## Overview
- **Material type:** `MultiMaterial`

## Parameters
### Mandatory
- `type`: Must be set to `"MultiMaterial"`
- `materials`: List of material names referenced as child blocks

Child material blocks must define their own parameters. Each child is any valid material model (e.g., `Isotropic`, `TransverseIsotropic`).
