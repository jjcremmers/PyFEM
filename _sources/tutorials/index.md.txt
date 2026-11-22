# Finite Element Tutorials

## Introduction

Most finite element textbooks focus on the mathematical derivation of the finite element method. Commercial software manuals, on the other hand, explain which buttons to press to perform an analysis.

These tutorials fill the gap between both worlds.

Using the Python finite element code accompanying this book, you will not only learn **how** to perform an analysis, but also **how a finite element program actually works internally**. Throughout the tutorials, you will inspect the source code, perform small modifications, and investigate the numerical algorithms behind commercial finite element software.

Each tutorial takes approximately **one hour** and consists of a combination of background theory, code exploration and practical exercises.

The tutorials assume that you have already studied the theoretical material presented in the book.

---

## Learning objectives

After completing the tutorial series you should be able to

- understand the complete workflow of a finite element program;
- interpret meshes, boundary conditions and solution vectors;
- perform reliable convergence studies;
- understand how finite element results are stored and visualised;
- understand the algorithms behind nonlinear analyses;
- understand the differences between implicit and explicit solvers;
- critically assess the quality and reliability of simulation results.

---

# Tutorial overview

## 1. Your First Finite Element Analysis

> Building and solving your first finite element model.

Topics

- Mesh generation
- Node numbering
- Element connectivity
- Assembly
- Solving the linear system

[Go to Tutorial 1](tutorial01.md)

---

## 2. Boundary Conditions and Loads

Learn how loads and constraints are actually implemented inside a finite element code.

Topics

- Essential and natural boundary conditions
- Concentrated and distributed loads
- Reaction forces
- Singular stiffness matrices

[Go to Tutorial 2](tutorial02.md)

---

## 3. Mesh Quality and Convergence

A simulation is only as good as its mesh.

Topics

- Mesh refinement
- Convergence studies
- Element quality
- Stress singularities

[Go to Tutorial 3](tutorial03.md)

---

## 4. Post-processing: Where Do the Results Come From?

Understand what commercial FE software is actually plotting.

Topics

- Nodal quantities
- Element quantities
- Integration points
- Stress recovery
- Averaging

[Go to Tutorial 4](tutorial04.md)

---

## 5. Geometrically Nonlinear Analysis

Learn how Newton-Raphson iterations solve nonlinear equilibrium.

Topics

- Residual forces
- Tangent stiffness
- Load stepping
- Convergence

[Go to Tutorial 5](tutorial05.md)

---

## 6. Material Nonlinearities

Introduce constitutive models and history-dependent materials.

Topics

- Plasticity
- History variables
- Return mapping
- Internal variables

[Go to Tutorial 6](tutorial06.md)

---

## 7. Contact and Constraints

How finite element programs detect and enforce contact.

Topics

- Contact detection
- Penalty methods
- Constraint equations
- Tied interfaces

[Go to Tutorial 7](tutorial07.md)

---

## 8. Dynamics and Explicit Analysis

Extend the solver to transient problems.

Topics

- Mass matrices
- Time integration
- Explicit solvers
- Stability

[Go to Tutorial 8](tutorial08.md)

---

## 9. Advanced Element Technology

Discover why different elements produce different answers.

Topics

- Reduced integration
- Locking
- Hourglassing
- Mixed formulations

[Go to Tutorial 9](tutorial09.md)

---

## 10. Verification, Validation and Debugging

How experienced analysts know whether they can trust a simulation.

Topics

- Patch tests
- Equilibrium
- Energy balance
- Error diagnosis

[Go to Tutorial 10](tutorial10.md)

---

## Suggested workflow

The tutorials are intended to be completed in order.

Each tutorial builds upon concepts introduced previously and extends the finite element code with new functionality.

Approximate time per tutorial:

**45–60 minutes**