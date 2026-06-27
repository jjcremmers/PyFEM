# Beam3D Examples

This directory contains small displacement-driven examples for the `Beam3D`
element.

## Files

- `cantilever_axial.dat` / `cantilever_axial.pro`: two-element cantilever with a prescribed axial tip displacement
- `cantilever_torsion.dat` / `cantilever_torsion.pro`: two-element cantilever with a prescribed tip twist
- `deployable_ring.dat` / `deployable_ring.pro`: simplified closed circular deployable-ring benchmark with three equally spaced tangential torsion actuators
- `vertical_member.dat` / `vertical_member.pro`: two-element vertical member with a prescribed axial shortening

## Run

```bash
pyfem cantilever_axial.pro
pyfem cantilever_torsion.pro
pyfem deployable_ring.pro
pyfem vertical_member.pro
```

These examples are intentionally simple and use prescribed displacements rather
than applied distributed loads so the element kinematics can be inspected
directly through the reported reactions. The deployable ring is the exception:
it uses three nodal moment pairs aligned with the local ring tangents so the
actuation is torsional rather than an arbitrary global rotation component.

`Beam3D` now accepts an optional `orientation = [x1, x2, x3]` property. This
vector defines the preferred local beam `x` direction; it must not be parallel
to the beam axis.

`Beam3D` also accepts an optional `sectionStiffness` property as a flat
36-entry row-major list defining a `6 x 6` section stiffness matrix. The
current element supports coupling within the axial/shear block and within the
curvature/torsion block, but not coupling between those two `3 x 3` blocks.

For the scalar-property path, `Beam3D` now also accepts the explicit section
stiffness entries `EA`, `GAy`, `GAz`, `EIy`, `EIz`, and `GJ`. When these are
not given, the element falls back to the older `E`, `G`, `A`, `Ix`, `Iy`, and
`J` properties.

For inertia, `Beam3D` now accepts `massPerLength`, `rotaryInertia1`,
`rotaryInertia2`, and `polarRotaryInertia`. The legacy fallback remains
`rho * A`, `rho * Ix`, `rho * Iy`, and `rho * J`. The aliases `rhoA`, `rhoIx`,
`rhoIy`, and `rhoJ` are also accepted.

## Reference

The deployable-ring example is inspired by:

Y. Goto, Y. Watanabe, T. Kasugai, and M. Obata,
“Elastic buckling phenomenon applicable to deployable rings,”
*International Journal of Solids and Structures* 29(7), 893-909, 1992.
DOI: `10.1016/0020-7683(92)90024-N`

The current PyFEM example is a reduced demonstrator of the same class of
problem. It is not intended as an exact reproduction of the full geometry,
boundary conditions, or loading path from the paper.
