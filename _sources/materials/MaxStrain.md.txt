# MaxStrain

`MaxStrain` is a composite-material failure criterion that compares the
normal and shear strains with allowable tensile, compressive, and shear strain
limits. Failure is indicated when the largest normalized strain reaches one.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Failure index

For the material axes, the criterion evaluates the fiber strain, transverse
strain, and in-plane shear strain independently. Tension and compression use
separate limits. The failure index is the maximum normalized value; an index
greater than or equal to 1 indicates failure.

## Input parameters

```text
failure =
{
  type    = "MaxStrain";
  eps1t   = 0.015;
  eps1c   = 0.012;
  eps2t   = 0.005;
  eps2c   = 0.010;
  gamma12u = 0.020;
};
```

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"MaxStrain"`. | String | Required |
| `eps1t`, `eps1c` | Tensile and compressive strain limits in direction 1. | Float | Required |
| `eps2t`, `eps2c` | Tensile and compressive strain limits in direction 2. | Float | Required |
| `gamma12u` | Allowable engineering shear strain. | Float | Required |

The strain components must be expressed in the material coordinate system.
