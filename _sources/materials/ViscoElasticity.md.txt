# ViscoElasticity

`ViscoElasticity` is a time-dependent linear viscoelastic material based on a
generalized Maxwell model. An equilibrium spring is combined with Maxwell
elements, each having a relaxation modulus and relaxation time. The material
therefore captures stress relaxation and creep in three-dimensional or
two-dimensional continuum analyses.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

```text
material =
{
  type        = "ViscoElasticity";
  E           = 1000.0;
  nu          = 0.3;
  Einf        = 100.0;
  nMaxwell    = 2;
  relaxTimes  = [1.0, 10.0];
  relaxModuli = [300.0, 600.0];
};
```

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"ViscoElasticity"`. | String | Required |
| `E` | Instantaneous Young's modulus. | Float | Required |
| `nu` | Poisson's ratio. | Float | Required |
| `Einf` | Long-term equilibrium modulus. | Float | Optional; default = `0.1*E` |
| `nMaxwell` | Number of Maxwell elements. | Integer | Optional; default = 3 |
| `relaxTimes` | Relaxation time for each Maxwell element. | List of floats | Optional; generated if omitted |
| `relaxModuli` | Relaxation modulus for each Maxwell element. | List of floats | Optional; equally distributed if omitted |

If the lists are omitted, default values are generated. Their lengths must
match `nMaxwell`. The model stores internal strains and requires time
increments from the solver.

## Examples

- [Creep test](../../examples/materials/viscoelasticity/creep_test.pro)
- [Strain-step test](../../examples/materials/viscoelasticity/strain_step.pro)
