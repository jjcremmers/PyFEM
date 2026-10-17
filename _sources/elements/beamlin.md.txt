# BeamLin

`BeamLin` is a linear two-node beam element for small-deformation structural
analyses. It uses the beam kinematics and material/section properties supplied
in the element group to compute axial, bending, and shear response.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Overview

- **Element type:** `BeamLin`
- **Application:** linear beam analysis
- **Optional output:** body-force loading through `bodyForce`

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"BeamLin"`. | String | Required |
| `material` | Beam material or section definition. | Block | Required by the analysis |
| `bodyForce` | Distributed body-force specification. | Float/list | Optional |

Verify the required nodal degrees of freedom and section convention against
the selected beam formulation.
