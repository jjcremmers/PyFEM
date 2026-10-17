# ThermoSmallStrainAxiSym

`ThermoSmallStrainAxiSym` is an axisymmetric coupled thermo-mechanical element
using a small-strain formulation. It is suitable for rotationally symmetric
thermal and mechanical loading.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"ThermoSmallStrainAxiSym"`. | String | Required |
| `rank` | Must be 2 for the axisymmetric formulation. | Integer | Required |
| `material` | Thermo-mechanical material block. | Block | Required |
