# FiniteStrainAxiSym

The `FiniteStrainAxiSym` element performs axisymmetric finite strain analysis for rotationally symmetric structures undergoing large deformations.

## Overview

The element uses axisymmetric finite-strain kinematics for bodies and loading
conditions with rotational symmetry. It supports geometric and material
nonlinearity, including large displacements and rotations.

## Parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Element formulation. | String | Required; must be `"FiniteStrainAxiSym"`. |
| `material` | Continuum material formulation for axisymmetric finite-strain analysis. | String | Required; the block must contain a material-model `type`. |

The selected material must provide a constitutive response compatible with the
finite-strain kinematics and axisymmetric stress state. See the [materials
overview](../materials/overview.md) for available formulations and required
parameters.
