# SmallStrainAxiSym

The `SmallStrainAxiSym` element performs axisymmetric small strain analysis for rotationally symmetric structures.

## Overview

The element uses axisymmetric small-strain kinematics for bodies and loading
conditions with rotational symmetry. It is suitable for linear analysis and
for constitutive models whose response remains within the small-strain
assumption.

## Parameters

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Element formulation. | String | Required; must be `"SmallStrainAxiSym"`. |
| `material` | Continuum material formulation for axisymmetric small-strain analysis. | Block | Required; the block must contain a material-model `type`. |

The selected material must provide a constitutive response compatible with the
small-strain kinematics and axisymmetric stress state. See the [materials
overview](../materials/overview.md) for available formulations and required
parameters.
