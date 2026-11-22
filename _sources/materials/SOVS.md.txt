# Skorohod–Olevsky Viscous Sintering

The `SkorohodOlevsky` model describes viscous densification and deformation of
ceramic powder during sintering. It combines volumetric flow driven by
sintering stress with deviatoric viscous deformation driven by applied stress.
The viscosity depends on temperature and relative density, allowing the model
to represent free sintering, pressure-assisted sintering, and creep during
densification.

```{warning}
This page was generated with the assistance of AI and must be checked
carefully against the implementation and the conventions of the element being
used.
```

## Input parameters

```text
material =
{
  type       = "SkorohodOlevsky";
  eta0       = 1.0e12;
  Q          = 500000.0;
  T          = 1600.0;
  rho0       = 0.6;
  sigma_sint = 1.0e6;
  n_vol      = 2.0;
  n_shear    = 1.0;
  R          = 8.314;
};
```

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"SkorohodOlevsky"`. | String | Required |
| `eta0` | Reference viscosity. | Float | Required |
| `Q` | Activation energy for viscous flow. | Float | Required |
| `T` | Sintering temperature. | Float | Required |
| `rho0` | Initial relative density in `(0, 1]`. | Float | Required |
| `sigma_sint` | Sintering stress. | Float | Required |
| `n_vol` | Volumetric viscosity exponent. | Float | Optional; default = 2.0 |
| `n_shear` | Shear viscosity exponent. | Float | Optional; default = 1.0 |
| `R` | Gas constant. | Float | Optional; default = 8.314 |

The model reports the stress components, current density, and volume strain.
It stores density, viscous strain, stress, and time history variables, so it
requires a time-dependent solver configuration.

## Examples

- [Free sintering](../../examples/materials/sintering/free_sintering.pro)
- [Pressure-assisted sintering](../../examples/materials/sintering/pressure_sintering.pro)
