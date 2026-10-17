# Tutorial: Using GMSH Input Files in PyFEM

This tutorial explains how to use GMSH-generated meshes as input for PyFEM, based on the examples in the [examples/gmsh](../../examples/gmsh) directory. We will walk through the structure of the input files and how to set up a simulation using GMSH meshes.

## 1. Overview of Files

A typical GMSH-based example in PyFEM consists of:

- `.geo` file: GMSH geometry definition
- `.msh` file: GMSH mesh file (exported from GMSH)
- `.dat` file: PyFEM input file referencing the mesh and specifying boundary conditions, loads, etc.
- `.pro` file: PyFEM project file (optional, for advanced setups)

## 2. Creating the Mesh with GMSH

Write a geometry file (e.g., [two_fibres.geo](../../examples/gmsh/two_fibres.geo)):

```geo
lc = 0.25;
Point(1)  = { 0.0, 0.0, 0.0, lc };
Point(2)  = { 10.0, 0.0, 0.0, lc };
// ... more points ...
Line(1) = {1,2};
// ... more lines, circles, surfaces ...
Physical Surface("Epoxy")  = {1};
Physical Surface("Fibre")  = {2};
Physical Surface("Fibre2") = {3};
```

Generate the mesh in GMSH:

```bash
gmsh two_fibres.geo -2 -format msh4 -o two_fibres.msh
```

## 3. PyFEM Input File Structure

The `.dat` file links the mesh and defines boundary conditions:

```ini
gmsh = "two_fibres.msh"

<NodeConstraints>
  u[Right] = 2.5;
  v[Right] = 0.0;
  u[Left]  = 0.0;
  v[Left]  = 0.0;
</NodeConstraints>

<ExternalForces>
</ExternalForces>
```

- `gmsh = "two_fibres.msh"` tells PyFEM to use the GMSH mesh (see [two_fibres.msh](../../examples/gmsh/two_fibres.msh)).
- `<NodeConstraints>` block applies Dirichlet boundary conditions to nodes in the physical groups `Right` and `Left` (as defined in the `.geo` file).
- `<ExternalForces>` block can be used to apply loads (empty in this example).

## 4. Running the Simulation

From the PyFEM root directory, run:

```bash
pyfem two_fibres.dat
```

## 5. Visualizing Results

PyFEM writes output files (e.g., `.pvd`, `.vtu`) for visualization in ParaView:

```bash
paraview two_fibres.pvd
```

## 6. More Examples

- [twist.geo](../../examples/gmsh/twist.geo), [twist.msh](../../examples/gmsh/twist.msh), [twist.dat](../../examples/gmsh/twist.dat) show a 3D example with extrusion and different boundary conditions.
- The process is the same: define geometry in `.geo`, mesh with GMSH, reference the mesh in the `.dat` file, and specify boundary conditions using physical group names.

## 7. Tips

- Always define physical groups in your `.geo` file for all regions where you want to apply boundary conditions or loads.
- Use GMSH's GUI to inspect and assign physical groups interactively.
- For advanced setups, see the `.pro` files (e.g., [twist.pro](../../examples/gmsh/twist.pro) or [two_fibres.pro](../../examples/gmsh/two_fibres.pro)) or the [PyFEM documentation](../../README.md).

---

For more details, see the [examples/gmsh](../../examples/gmsh) directory and the PyFEM documentation.
