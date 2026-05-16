# Reissner-Mindlin Shell Example

This example drives a curved cylindrical shell strip with the `ReissnerMindlinShell` element in a nonlinear analysis.

## Files

- `curved_cantilever.dat`: curved cylindrical strip with prescribed end displacement
- `curved_cantilever.pro`: solver and output settings for the curved cantilever
- `curved_cantilever_composite.pro`: same curved cantilever using a layered composite laminate
- `pinched_hemisphere.dat`: quarter pinched hemisphere with an 18 degree hole
- `pinched_hemisphere.pro`: solver and output settings for the hemisphere benchmark
- `pinched_hemisphere_fine.dat`: finer quarter hemisphere mesh ported from the `dawn` shell benchmark
- `pinched_hemisphere_fine.pro`: solver and output settings for the fine hemisphere benchmark

## Model

- Geometry: cylindrical strip with radius about `2.0`, length `4.0`, width `1.5592`
- Mesh: `4` curved Quad4 shell elements
- Boundary conditions: left edge clamped, right edge prescribed downward displacement in `w`
- Solver: `NonlinearSolver` with fixed load steps

## Run

```bash
pyfem curved_cantilever.pro
pyfem curved_cantilever_composite.pro
pyfem pinched_hemisphere.pro
pyfem pinched_hemisphere_fine.pro
```

The graph output reports the prescribed displacement at node `5` and the corresponding internal reaction in `w`.
For the hemisphere case, the graph output reports the loaded-point displacement and reaction in `u` at node `72`.
For the fine hemisphere port, the graph output reports the two symmetry-edge pinch displacements
at nodes `1` and `17` together with the load factor `lam`.

## Composite Laminate Input

The shell element accepts the same laminate definition pattern used by `Plate`:

```text
Shell =
{
  type = "ReissnerMindlinShell";

  materials = [ "UD" ];
  layers    = [ "ply0_bot" , "ply90_bot" , "ply90_top" , "ply0_top" ];

  UD =
  {
    E1   = 1.35e5;
    E2   = 1.0e4;
    nu12 = 0.3;
    G12  = 5.0e3;
    G13  = 4.0e3;
    G23  = 3.8e3;
    rho  = 1.6e-9;
  };
};
```

Each layer block sets `material`, `theta`, and `thickness`. The shell integrates each ply separately through the thickness using the laminate data from `pyfem.elements.Composite.Laminate`.
