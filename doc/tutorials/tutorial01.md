# Tutorial 1 – Your First Finite Element Analysis

## Estimated time

60 minutes

---

## Learning objectives

After completing this tutorial you should be able to

- explain the complete finite element workflow;
- describe the purpose of nodes and elements;
- understand the global stiffness matrix;
- solve a small linear finite element problem.

---

## Background

This tutorial introduces the basic workflow of every finite element program.

Although commercial software performs these steps automatically, every finite element analysis follows exactly the same sequence.

1. Generate a mesh
2. Create elements
3. Assemble the global stiffness matrix
4. Apply boundary conditions
5. Solve the linear system
6. Post-process the results

---

## Inside the code

Files to inspect

- mesh.py
- assembly.py
- solver.py

Questions

- Where is the global stiffness matrix created?
- Where are the element matrices assembled?
- How many degrees of freedom does the model contain?

---

## Exercise 1

Run the example problem.

Inspect

- number of nodes
- number of elements
- size of the stiffness matrix

---

## Exercise 2

Double the number of elements.

How does the solution change?

---

## Things to investigate

- What happens if you remove one element?
- What happens if two nodes occupy the same location?
- Does node numbering influence the solution?

---

## Common pitfalls

- Confusing nodes with elements.
- Assuming the mesh stores stresses.
- Forgetting that each node contains multiple degrees of freedom.

---

## Summary

You now understand the complete workflow of a linear finite element analysis.

Next tutorial:

Boundary conditions and loads.