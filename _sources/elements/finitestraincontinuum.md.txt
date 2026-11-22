# FiniteStrainContinuum

The `FiniteStrainContinuum` element is a displacement-based finite element for 2D and 3D large deformation analysis in continuum mechanics. This element formulation accounts for geometric nonlinearity through finite strain kinematics, making it suitable for problems involving large displacements, large rotations, and moderate to large strains.

## Overview

The element supports total-Lagrangian (TL) and updated-Lagrangian (UL)
formulations in two and three dimensions. It accounts for geometric and
material nonlinearity and is intended for large displacements, rotations, and
strains.

## Parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Element formulation. | String | Required; must be `"FiniteStrainContinuum"`. |
| `material` | Continuum material formulation used for the constitutive response. | Block | Required; the block must contain a material-model `type`. |
| `method` | Finite-strain kinematic formulation. | String | Optional; `"TL"` selects Total Lagrangian and `"UL"` selects Updated Lagrangian. Defaults to `"TL"`. |

The material model and its parameters depend on the finite-strain
constitutive formulation. Consult the [materials overview](../materials/overview.md)
and the relevant material page for required properties, assumptions, and
valid strain measures.
