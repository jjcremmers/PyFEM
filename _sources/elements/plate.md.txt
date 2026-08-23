# Plate

The plate element implements the Kirchhoff-Love finite element for flat structures in the x-y plane as 3, 4, 6, and 8 node elements with 5 degrees of freedom per node (u, v, w, rx, ry).

## Parameters
### Mandatory
- `type`: Must be set to `"Plate"`
- `material`: Single-material block defining constitutive behavior (e.g., `E`, `nu`, `rho`)
- `thickness`: Plate thickness

### Optional
(See documentation for details)
