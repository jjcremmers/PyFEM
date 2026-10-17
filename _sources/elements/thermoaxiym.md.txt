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
| `type` | Must be `"ThermoAxiSym"`. | String | Required |
| `rank` | Must match the axisymmetric formulation. | Integer | Required |
| `material` | Thermo-mechanical material block. | Block | Required |

Use this element only when the geometry, material orientation, and loading are
axisymmetric.
