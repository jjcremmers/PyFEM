# PhaseField

`PhaseField` is a coupled element for phase-field fracture analyses. It
combines the displacement field with a phase or damage field and uses a
constitutive material to evaluate the mechanical response. The phase field
regularizes crack development over a finite length scale.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"PhaseField"`. | String | Required |
| `rank` | Spatial dimension of the formulation. | Integer | Required; typically 2 or 3 |
| `material` | Mechanical constitutive material block. | Block | Required |

The phase-field parameters are formulation- and material-dependent. Check the
selected material and example input file carefully.

## Examples

- [Phase-field patch test](../../examples/elements/phasefield/PhaseFieldPatch4.pro)
- [Phase-field dogbone](../../examples/elements/phasefield/dogbone/Dogbone9.pro)
