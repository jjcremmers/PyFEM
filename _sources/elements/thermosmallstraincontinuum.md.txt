# ThermoSmallStrainContinuum

`ThermoSmallStrainContinuum` is a coupled thermo-mechanical continuum element
for small-strain analyses. It solves the displacement and temperature fields
using a small-strain mechanical formulation.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"ThermoSmallStrainContinuum"`. | String | Required |
| `rank` | Spatial dimension. | Integer | Required |
| `material` | Thermo-mechanical material block. | Block | Required |

Thermal expansion and conductivity data must be provided consistently by the
material and analysis setup.
