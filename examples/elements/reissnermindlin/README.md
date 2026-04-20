# Reissner-Mindlin Shell Example

This example drives a curved cylindrical shell strip with the `ReissnerMindlinShell` element in a nonlinear analysis.

## Files

- `curved_cantilever.dat`: curved cylindrical strip with prescribed end displacement
- `curved_cantilever.pro`: solver and output settings for the curved cantilever
- `pinched_hemisphere.dat`: quarter pinched hemisphere with an 18 degree hole
- `pinched_hemisphere.pro`: solver and output settings for the hemisphere benchmark

## Model

- Geometry: cylindrical strip with radius about `2.0`, length `4.0`, width `1.5592`
- Mesh: `4` curved Quad4 shell elements
- Boundary conditions: left edge clamped, right edge prescribed downward displacement in `w`
- Solver: `NonlinearSolver` with fixed load steps

## Run

```bash
pyfem curved_cantilever.pro
pyfem pinched_hemisphere.pro
```

The graph output reports the prescribed displacement at node `5` and the corresponding internal reaction in `w`.
For the hemisphere case, the graph output reports the loaded-point displacement and reaction in `u` at node `72`.
