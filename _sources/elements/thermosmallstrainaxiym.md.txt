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
| `type` | Element formulation. | String | Required; must be `"ThermoSmallStrainAxiSym"`. |
| `rank` | Must be 2 for the axisymmetric formulation. | Integer | Required |
| `material` | Continuum thermo-mechanical material formulation. | Block | Required; the block must contain a material-model `type`. |

See the [materials overview](../materials/overview.md) for the available
thermo-mechanical material formulations and their required parameters.
