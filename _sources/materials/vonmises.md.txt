# VonMises

`VonMises` failure index computes the von Mises stress-based criterion compared to an allowable stress `smax`.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Overview
- **Material type:** `VonMises` (failure criterion)

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"VonMises"`. | String | Required |
| `smax` | Allowable von Mises stress. | Float | Required |

## Examples
- Post-processing failure checks in continuum elements
