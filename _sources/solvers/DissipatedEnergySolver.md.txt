# DissipatedEnergySolver

## Background

The `DissipatedEnergySolver` is a path-following solver for nonlinear
problems with material failure, softening, delamination, or other phenomena
that can make conventional load control difficult. It follows the equilibrium
path through limit points by controlling the energy dissipated during an
increment.

The method is particularly useful when the load–displacement response contains
snap-back or unstable branches. It is described in:

> E. Börjesson, J.J.C. Remmers, and M. Fagerström, “A generalised
> path-following solver for robust analysis of material failure,”
> *Computational Mechanics* (2023).
>
> [doi:10.1007/s00466-022-02175-w](https://doi.org/10.1007/s00466-022-02175-w)

The solver requires elements or material models that provide a dissipation
contribution. During assembly, PyFEM collects the element dissipation vector
and the total dissipated energy. Interface elements with cohesive material
models and continuum elements with inelastic material models are typical
applications.

## Implementation

The implementation is in `pyfem/solvers/DissipatedEnergySolver.py`. Each load
cycle uses Newton–Raphson iterations. The solver assembles:

- the tangent stiffness and internal force;
- the external reference load vector;
- the dissipation gradient and the accumulated dissipated energy.

The solver starts with an arc-length-controlled predictor–corrector step. Once
the dissipated energy in a converged step exceeds `switchEnergy`, it switches
to dissipated-energy control. In the default `Local` formulation, the energy
constraint is based directly on the assembled dissipation:

```{math}
g = \mathcal{D} - \Delta\tau,
```

where $\mathcal{D}$ is the current dissipated energy and $\Delta\tau$ is the
target energy increment. The corresponding dissipation gradient is used in the
Newton correction.

The alternative `Classic` formulation evaluates the energy constraint from
the load factor and displacement increments. The solver adapts the target
energy increment using the number of Newton iterations. The adaptation aims
to approach `optiter` iterations and is limited by `maxdTau`.

After convergence, the solver stores the converged increment, updates the
target energy increment, commits element history, and stops when either
`maxLam` or `maxCycle` is reached.

:::{warning}
`switchEnergy` is required by the implementation. A model must also contain
elements or materials that implement `getDissipation`; otherwise the
dissipation contribution remains zero and energy control cannot represent the
intended failure process.
:::

## Input parameters

The solver is selected with:

```text
solver =
{
  type = "DissipatedEnergySolver";
};
```

The available solver parameters are:

| Parameter | Description | Type | Remarks |
| --- | --- | --- | --- |
| `type` | Must be `"DissipatedEnergySolver"`. | String | — |
| `tol` | Relative residual tolerance for Newton–Raphson convergence. | Float | Optional, default = $10^{-4}$ |
| `optiter` | Target number of Newton iterations used for step-size adaptation. | Int | Optional, default = 5 |
| `iterMax` | Maximum Newton iterations. | Int | Optional, default = 10 |
| `maxdTau` | Upper bound for the target dissipated-energy increment. | Float | Optional, default = $10^{20}$ |
| `maxLam` | Maximum load factor. | Float | Optional, default = $10^{20}$ |
| `lam` | Initial load factor. |Float | Optional, default = 1.0 |
| `disstype` | Dissipation constraint formulation: `"Local"` or `"Classic"`. | String | Optional, default = '"Local"` |
| `switchEnergy` | Dissipated energy at which the solver switches from arc-length to energy control. | Float | — |
| `maxCycle` | Maximum number of load cycles. | Int | Optional, default = 1000 |

The solver starts internally in `arclength-controlled` mode. The current
implementation does not expose the initial control method as an input option;
`switchEnergy` determines when energy control begins.

### Dissipation support

Elements contribute dissipation through a `getDissipation` method. The
assembled contribution is based on the element's current state and material
history. For example, the `Interface` element and continuum elements with
inelastic constitutive models can contribute to the solver's energy measure.

## Examples

The reference examples are:

- [Delamination buckling with 100 elements](../../examples/solver/dissipatedEnergySolver/delam_buckling100.pro)
  uses a finite-strain continuum together with interface elements.
- [Delamination buckling with 200 elements](../../examples/solver/dissipatedEnergySolver/delam_buckling200.pro)
  uses the same model with a finer mesh.
- [Interface peel test](../../examples/elements/interface/PeelTest.pro)
  demonstrates local dissipation control for a cohesive interface.
- [Peel test from Chapter 13](../../examples/ch13/PeelTest60.pro)
  combines the solver with graph, contour, and HDF5 output.

A typical solver block is:

```text
solver =
{
  type         = "DissipatedEnergySolver";
  maxCycle     = 60;
  tol          = 1.0e-3;
  maxLam       = 50;
  lam          = 1.0;
  disstype     = "Local";
  switchEnergy = 1.0e-3;
  maxdTau      = 0.05;
};
```

The detailed tutorial for the delamination buckling case is available in
[The dissipated energy path-following solver](../tutorials/dissnrg_tutorial.md).
