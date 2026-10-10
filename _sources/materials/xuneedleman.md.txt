# Xu–Needleman

`XuNeedleman` is a cohesive-zone material for interface elements. It relates
the displacement jump across an interface to normal and tangential tractions,
and is suitable for modelling debonding, delamination, and peeling.

The model is characterized by two physical parameters:

- `Tult`: the ultimate traction;
- `Gc`: the fracture energy.

The material supports both two-dimensional interfaces, with one normal and one
tangential component, and three-dimensional interfaces, with one normal and
two tangential components.

```{warning}
`XuNeedleman` is an interface material. It must be assigned to an interface
element, such as `Interface`; it is not a bulk continuum material.
```

## Implementation

The implementation is in `pyfem/materials/XuNeedleman.py`. For a two-dimensional
interface, the deformation vector is interpreted as

```{math}
\boldsymbol{v} =
\begin{bmatrix}
v_n \\
v_t
\end{bmatrix},
```

and for a three-dimensional interface as

```{math}
\boldsymbol{v} =
\begin{bmatrix}
v_n \\
v_{t1} \\
v_{t2}
\end{bmatrix}.
```

The characteristic normal and tangential separation scales are calculated
from `Gc` and `Tult`:

```{math}
v_{n,\max} = \frac{G_c}{e\,T_{\mathrm{ult}}},
\qquad
v_{t,\max} = \frac{G_c}{1.16580058\,T_{\mathrm{ult}}}.
```

The implementation uses fixed shape parameters $r = 0$ and $q = 1$. The
traction vector and consistent tangent are evaluated from the current
separation. In three dimensions, the tangential traction is resolved along
the direction of the tangential separation.

The cohesive potential is evaluated as

```{math}
\Phi(\boldsymbol{v}) = G_c + G_c e^{-v_n/v_{n,\max}}
\left[
\frac{(1-r+v_n/v_{n,\max})(1-q)}{r-1}
 - \left(q + \frac{r-q}{r-1}\frac{v_n}{v_{n,\max}}\right)
 e^{-(v_t/v_{t,\max})^2}
\right],
```

for nonnegative normal separation. The material tracks the dissipated energy
as a history variable. The incremental dissipation passed to the solver is
based on the difference between the cohesive potential and the current
traction–separation work:

```{math}
g = \Phi(\boldsymbol{v})
 - \frac{1}{2}\boldsymbol{v}\!\cdot\!\boldsymbol{t}
 - \Phi_{\mathrm{history}}.
```

The material stores the history only after the computed dissipation becomes
nonnegative. It is committed by the analysis framework after a converged load
step.

## Input parameters

The material is defined inside an interface element. A minimal configuration
is:

```text
InterfaceElem =
{
  type = "Interface";

  material =
  {
    type = "XuNeedleman";
    Tult = 0.5;
    Gc   = 0.1;
  };
};
```

The material parameters are:

| Parameter | Description | TYpe | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"XuNeedleman"`. | String | — |
| `Tult` | Ultimate cohesive traction. | Float | — |
| `Gc` | Fracture energy. | Float | — |

The following internal model constants are fixed by the implementation and
cannot currently be set through the input file:

| Parameter | Description | Value |
| --- | --- | --- |
| $r$ | Normal/tangential potential parameter. | $0$ |
| $q$ | Normal/tangential potential parameter. | $1$ |
| $v_{n,\max}$ | Characteristic normal separation. | $G_c/(eT_{\mathrm{ult}})$ |
| $v_{t,\max}$ | Characteristic tangential separation. | $G_c/(1.16580058T_{\mathrm{ult}})$ |

The material produces the following output labels:

| Interface rank | Output labels |
| --- | --- |
| 2D | `Tn`, `Ts` |
| 3D | `Tn`, `Ts1`, `Ts2` |

Units must be consistent. For example, if traction has units of force per
area and separation has units of length, `Gc` must have units of force per
length.

## Examples

The reference examples are:

- [Interface peel test](../../examples/elements/interface/PeelTest.pro)
  demonstrates a two-dimensional cohesive interface with `Tult = 0.5` and
  `Gc = 0.1`.
- [Chapter 13 peel test](../../examples/ch13/PeelTest60.pro) combines the
  material with a continuum domain, dissipated-energy path following, graph
  output, contour output, and HDF5 output.
- [Delamination buckling, 100 elements](../../examples/solver/dissipatedEnergySolver/delam_buckling100.pro)
  and [200 elements](../../examples/solver/dissipatedEnergySolver/delam_buckling200.pro)
  use `XuNeedleman` interface layers in a larger failure analysis.

The cohesive material is commonly paired with
`DissipatedEnergySolver`, which uses the material's dissipation contribution
to follow unstable failure paths. See the
[dissipated-energy solver documentation](../solvers/DissipatedEnergySolver.md)
for solver parameters and path-following details.
