# Contact Mechanics Model

## Overview
The contact model in PyFEM implements a penalty-based contact algorithm for simulating mechanical interactions between deformable bodies and rigid obstacles. The model supports both 2D (disc) and 3D (sphere) contact geometries with a simple penalty method that enforces non-penetration constraints.

The contact model is implemented as a special model type that can be added to any structural analysis to impose contact constraints between the mesh and a moving rigid obstacle.

## Theory
### Penalty Method
The contact algorithm uses the penalty method to enforce contact constraints. When a node penetrates a rigid obstacle, a penalty force is applied:

$$
\mathbf{f}_c = k_p \, g \, \mathbf{n}
$$
where:
- $k_p$ is the penalty stiffness parameter
- $g$ is the penetration depth (overlap)
- $\mathbf{n}$ is the outward normal direction from the obstacle surface

### Contact Detection
For a rigid obstacle (disc in 2D or sphere in 3D) with center $\mathbf{c}$ and radius $R$, contact is detected when:

$$
g = R - \|\mathbf{x} - \mathbf{c}\| > 0
$$
where $\mathbf{x}$ is the current position of a node. The penetration depth $g$ measures how far the node has penetrated into the obstacle.

The outward normal from the obstacle surface is:

(See documentation for further details)
