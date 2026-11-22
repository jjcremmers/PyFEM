# ThermoSurface

`ThermoSurface` is a thermal surface element for applying or modeling thermal
effects on a boundary. It is used with the thermal field and surrounding
thermo-mechanical formulation to represent surface heat transfer or thermal
boundary contributions.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"ThermoSurface"`. | String | Required |
| `material` | Thermal material or boundary-property block. | Block | Required by the analysis |

Check the element's boundary-node ordering and thermal boundary convention
against the selected analysis setup.
