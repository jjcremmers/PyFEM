# ViscoPlasticity

`ViscoPlasticity` is a rate-dependent plasticity model based on Perzyna
overstress theory. It combines isotropic elasticity with von Mises yielding,
viscoplastic flow above the yield stress, and optional linear hardening.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

```text
material =
{
  type   = "ViscoPlasticity";
  E      = 200000.0;
  nu     = 0.3;
  syield = 250.0;
  gamma  = 0.001;
  n      = 1.0;
  hard   = 1000.0;
};
```

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"ViscoPlasticity"`. | String | Required |
| `E` | Young's modulus. | Float | Required |
| `nu` | Poisson's ratio. | Float | Required |
| `syield` | Initial von Mises yield stress. | Float | Required |
| `gamma` | Fluidity parameter controlling rate sensitivity. | Float | Required |
| `n` | Overstress exponent. | Float | Optional; default = 1.0 |
| `hard` | Linear hardening modulus. | Float | Optional; default = 0.0 |

The model integrates the viscoplastic strain using the solver time increment
and stores elastic, plastic, and equivalent-plastic strain history variables.

## Examples

- [Bar tension with displacement control](../../examples/materials/viscoplasticity/bar_tension_dispcontrol.pro)
- [Bar tension with force control](../../examples/materials/viscoplasticity/bar_tension_forcecontrol.pro)
