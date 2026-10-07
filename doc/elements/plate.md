# Plate

`Plate` is a flat, small-displacement plate element for meshes in the global
`x-y` plane. Each node has five degrees of freedom: `u`, `v`, `w`, `rx`, and
`ry`. The element supports 4-node elements. Other surface elements (3-, 6- and 8-node 
elements) can be used, subject to the
shape functions used by the mesh. However, this is not tested.

The element accepts either a single material or a laminate. Layer angles are
given in degrees. Layers are integrated about the laminate mid-plane, so their
ordered thicknesses define the total plate thickness.

## Implementation

The implementation is in `pyfem/elements/Plate.py` and uses `Laminate` from
`pyfem/elements/Composite.py`. At each integration point it computes membrane
strains from `u` and `v`, curvatures from `rx` and `ry`, and transverse shear
strains from `w`, `rx`, and `ry`.

The membrane resultants and bending moments use the laminate constitutive
matrices:

```text
[N]   [A B] [eps0]
[M] = [B D] [kappa]
```

`A` is the extensional stiffness, `B` is the membrane-bending coupling
stiffness, and `D` is the bending stiffness. Transverse shear uses the
laminate shear matrix and a default shear-correction factor of `5/6`.

The element also provides a consistent mass matrix. Mass and rotary inertia
are calculated from the layer densities and layer positions; therefore `rho`
is required for dynamic analyses.

Stress output is evaluated at the bottom and top surfaces with these labels:

```text
s11bot, s22bot, s12bot, s11top, s22top, s12top
```

## User options (data)

The element block must contain `type = "Plate"` and one of the following
material definitions.

### Single material

Use `material` and `thickness` for a single-layer plate:

```text
PlateElem =
{
  type = "Plate";
  material =
  {
    type = "PlaneStress";  # optional material type tag
    E       = 1.0e6;
    nu      = 0.25;
    rho     = 1.0e3;
  };
  thickness = 0.1;
};
```

The material properties consumed by the plate are:

- `E`: Young's modulus. A scalar gives `E1 = E2`; a two-entry list gives
  `[E1, E2]`.
- `nu` or `nu12`: Poisson's ratio in the material axes.
- `G12`: in-plane shear modulus. If omitted, it is computed as
  `E1 / (2 * (1 + nu12))`.
- `G13` and `G23`: optional transverse shear moduli. If omitted, `G12` is
  used for both directions.
- `rho`: density.
- `thickness`: layer/plate thickness.

For this form, the layer angle defaults to `0` degrees.

### Laminate

For multiple materials or layers, define the material names in `materials`
and the layer names in `layers`:

```text
PlateElem =
{
  type      = "Plate";
  materials = ["UD"];
  layers    = ["l0", "l90", "l0"];

  UD =
  {
    E1   = 1.0e6;
    E2   = 5.0e5;
    nu12 = 0.25;
    G12  = 4.0e5;
    rho  = 1.0e3;
  };

  l0 = { material = "UD"; theta = 0.0;  thickness = 0.05; };
  l90 = { material = "UD"; theta = 90.0; thickness = 0.05; };
};
```

Each layer requires `material` and `thickness`; `theta` is optional and
defaults to `0` degrees. The order in `layers` is the order through the
thickness, from bottom to top. Total thickness is the sum of the layer
thicknesses.

The optional element property `shearCorrection` overrides the default `5/6`:

```text
shearCorrection = 0.8333333333;
```

Some older examples contain a `stack` entry. The current implementation does
not read `stack`; use the ordered `layers` list to define the laminate.

## Examples

- [Single-material cantilever](../../examples/elements/plate/plate_cantilever01.pro)
  demonstrates the basic `material` and `thickness` form.
- [Layered cantilever](../../examples/elements/plate/plate_cantilever02.pro)
  demonstrates a three-layer orthotropic laminate.
- [Plate tests](../../examples/plate/plate_test_01.pro) and the other
  `examples/plate/plate_test_*.pro` files show alternative material and layer
  definitions.
- [Dynamic plate](../../examples/plate/platedyn.pro) demonstrates eigenvalue
  analysis and the use of density.

A minimal static model is:

```text
input = "plate_cantilever01.dat";

PlateElem =
{
  type = "Plate";
  material = { E = 1.0e6; nu = 0.25; rho = 1.0e3; };
  thickness = 0.1;
};

solver = { type = "LinearSolver"; };
```
