# ThermoAxiSym

`ThermoAxiSym` is an axisymmetric thermo-mechanical continuum element. It is
intended for bodies and loading conditions with rotational symmetry and
couples axisymmetric displacement and temperature fields.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Element formulation. | String | Required; must be `"ThermoAxiSym"`. |
| `rank` | Must match the axisymmetric formulation. | Integer | Required |
| `material` | Continuum thermo-mechanical material formulation. | Block | Required; the block must contain a material-model `type`. |

Use this element only when the geometry, material orientation, and loading are
axisymmetric. See the [materials overview](../materials/overview.md) for the
available material formulations and their required parameters.
