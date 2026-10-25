# ThermoContinuum

`ThermoContinuum` is a continuum element for coupled thermo-mechanical
analyses. It combines displacement and temperature fields and passes the
thermal state to the selected material model.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"ThermoContinuum"`. | String | Required |
| `rank` | Spatial dimension. | Integer | Required |
| `material` | Thermo-mechanical material block. | Block | Required |

Thermal properties and field definitions are supplied through the selected
material and analysis configuration.
