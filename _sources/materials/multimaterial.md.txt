# MultiMaterial

The `MultiMaterial` model groups multiple material models for layered or zone-based analyses. It delegates `getStress()` to the active sub-material specified by the element (e.g., via `iMat`).

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Overview
- **Material type:** `MultiMaterial`

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"MultiMaterial"`. | String | Required |
| `materials` | List of material names referenced as child blocks. | List | Required |

Child material blocks must define their own parameters. Each child is any valid material model (e.g., `Isotropic`, `TransverseIsotropic`).
