# SolidLikeShell

`SolidLikeShell` is a shell-like continuum element that uses a solid-style
through-thickness formulation. It is intended for structures where shell
behavior is required while retaining a continuum representation of the
thickness.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"SolidLikeShell"`. | String | Required |
| `material` | Constitutive material block. | Block | Required |
| `thickness` | Shell or layer thickness. | Float | Required by the formulation |

Use a compatible continuum material and verify the element's nodal
kinematics before applying it to thin structures.
