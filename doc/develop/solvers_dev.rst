=======================
Developing Solver Modules
=======================

This guide explains how to implement new solution algorithms in PyFEM. Solvers
orchestrate the finite element analysis process, managing load steps, iteration
procedures, and convergence criteria to solve the discrete equilibrium equations.

Overview
--------

Solvers in PyFEM are responsible for:

- Managing the solution process (load steps, time integration)
- Implementing iteration schemes (Newton-Raphson, arc-length, explicit)
- Checking convergence criteria
- Calling assembly routines for system matrices
- Invoking output modules
- Managing solution status and history

All solver classes inherit from ``BaseModule`` and implement a ``run`` method
that performs one solution step (cycle).

Solver Class Structure
----------------------

Base Class
~~~~~~~~~~

All solvers inherit from ``BaseModule`` located in ``pyfem/util/BaseModule.py``.
This provides:

- Property management from input files
- Integration with the output system
- Logging capabilities

Required Methods
~~~~~~~~~~~~~~~~

A solver must implement:

.. code-block:: python

   class MySolver(BaseModule):
   
       def __init__(self, props, globdat):
           """Initialize solver with properties and global data."""
           BaseModule.__init__(self, props)
           # Initialize solver parameters
           
       def run(self, props, globdat):
           """Execute one solution step.
           
           Args:
               props: Properties tree (input file configuration)
               globdat: Global data container
               
           Returns:
               bool: True if analysis should continue
           """
           # Implementation here
           return continue_analysis

Implementation Examples
-----------------------

Example 1: Nonlinear Solver (Newton-Raphson)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows the pseudo code for a Newton-Raphson solver for nonlinear
static analysis, following Chapter 2 of the book *"Non-Linear Finite Element 
Analysis of Solids and Structures"* by de Borst et al.

The Newton-Raphson iteration solves (equation 2.21):

.. math::

   \mathbf{K}_t^{(i)} \Delta \mathbf{a}^{(i)} = \mathbf{f}_{ext} - \mathbf{f}_{int}^{(i)}

**Pseudo Code:**

.. code-block:: text

   class NonlinearSolver(BaseModule):
   
       def __init__(props, globdat):
           # Read parameters from input file
           tol = 1.0e-3
           iterMax = 10
           maxCycle = user_defined
           dtime = 1.0
           
           # Initialize load factor
           globdat.lam = 0.0
           
       def run(props, globdat):
           """Execute one load step using Newton-Raphson iteration.
           
           Algorithm from Box 2.1 (page 36):
           1. Apply load increment
           2. Iterate until convergence:
              a. Assemble tangent stiffness and internal forces
              b. Compute residual
              c. Solve for displacement increment
              d. Update displacement
              e. Check convergence
           """
           
           # Initialize displacement increment
           Da = zeros(dofCount)
           
           # Update load factor
           globdat.lam = loadFunction(cycle * dtime)
           
           # Assemble external force vector
           fext = assembleExternalForce(props, globdat)
           
           # Newton-Raphson iteration loop
           error = 1.0
           iiter = 0
           
           while error > tol:
               iiter += 1
               if iiter > iterMax:
                   raise RuntimeError("Not converged")
               
               # Assemble tangent stiffness and internal forces
               K, fint = assembleTangentStiffness(props, globdat)
               
               # Compute residual: r = f_ext - f_int
               residual = fext - fint
               
               # Apply boundary conditions
               K_constrained, residual_constrained = applyConstraints(K, residual)
               
               # Solve linear system: K * da = residual
               da = solve(K_constrained, residual_constrained)
               
               # Update displacement
               Da += da
               state += da
               
               # Compute convergence error
               error = norm(residual) / norm(fext)
           
           # Check stopping criteria
           if cycle >= maxCycle or lam >= maxLam:
               return False  # Stop analysis
           
           return True  # Continue to next step

**Full Implementation:**

See the complete implementation in:
`pyfem/solvers/NonlinearSolver.py <https://github.com/jjcremmers/PyFEM/blob/master/pyfem/solvers/NonlinearSolver.py>`_

Example 2: Arc-Length Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For problems with snap-back or snap-through behavior, implement arc-length
control following Chapter 2.5 of the book.

The Riks-Wempner arc-length method (equation 2.71):

.. math::

   \Delta \mathbf{a}^T \Delta \mathbf{a} + \psi^2 \Delta \lambda^2 = \Delta l^2

**Pseudo Code:**

.. code-block:: text

   class ArcLengthSolver(BaseModule):
   
       def __init__(props, globdat):
           # Read parameters
           arcLength = 1.0
           tol = 1.0e-3
           iterMax = 10
           psi = 1.0  # Load-displacement ratio
           
           # Initialize load factor
           globdat.lam = 0.0
   
       def run(props, globdat):
           """Execute one step with arc-length control.
           
           Algorithm from Section 2.5:
           1. Predictor: compute displacement direction
           2. Scale by arc-length constraint
           3. Corrector iterations with constraint
           """
           
           # Assemble reference load vector (unscaled)
           fref = assembleExternalForce(props, globdat)
           
           # PREDICTOR STEP
           K, fint = assembleTangentStiffness(props, globdat)
           
           # Solve for displacement direction
           da1 = solve(K, fref)
           
           # Compute initial load increment from arc-length (eq. 2.72)
           Dlam = arcLength / sqrt(da1^T * da1 + psi^2)
           
           # Apply predictor
           Da = Dlam * da1
           state += Da
           lam += Dlam
           
           # CORRECTOR ITERATIONS
           error = 1.0
           iiter = 0
           
           while error > tol:
               iiter += 1
               if iiter > iterMax:
                   raise RuntimeError("Not converged")
               
               # Assemble at current state
               K, fint = assembleTangentStiffness(props, globdat)
               
               # Compute residual
               fext = lam * fref
               residual = fext - fint
               
               # Solve for two correction vectors
               da_r = solve(K, residual)    # Residual direction
               da_f = solve(K, fref)        # Load direction
               
               # Arc-length constraint (equation 2.71)
               # Solve quadratic equation for load multiplier
               A = da_f^T * da_f + psi^2
               B = 2.0 * (Da^T * da_f + Da^T * da_r)
               C = Da^T * Da - arcLength^2
               
               dlam = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
               
               # Compute displacement correction
               da = da_r + dlam * da_f
               
               # Update
               Da += da
               state += da
               lam += dlam
               
               # Check convergence
               error = norm(residual) / norm(fext)
           
           # Adapt arc-length for next step based on iteration count
           if iiter < 4:
               arcLength *= 1.5
           elif iiter > 7:
               arcLength *= 0.5
           
           return True

**Full Implementation:**

See the complete implementation in:
`pyfem/solvers/RiksSolver.py <https://github.com/jjcremmers/PyFEM/blob/master/pyfem/solvers/RiksSolver.py>`_

Example 3: Explicit Dynamics Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For dynamic problems, implement explicit time integration following Chapter 9.

The central difference method (equation 9.19) integrates:

.. math::

   \mathbf{M} \mathbf{a}_{n+1} = \mathbf{f}_{ext} - \mathbf{f}_{int}(\mathbf{u}_n)

**Pseudo Code:**

.. code-block:: text

   class ExplicitSolver(BaseModule):
   
       def __init__(props, globdat):
           # Read parameters
           dt = user_defined  # Must satisfy CFL condition
           endTime = 1.0
           rho = material_density
           
           # Initialize state vectors
           globdat.velo = zeros(dofCount)
           globdat.acce = zeros(dofCount)
           
           # Assemble mass matrix (constant for explicit)
           M = assembleMassMatrix(props, globdat)
           
           # Compute lumped mass (diagonal approximation)
           Mlumped = sum(M, axis=1)
           
           globdat.time = 0.0
   
       def run(props, globdat):
           """Execute one explicit time step.
           
           Central difference scheme (Box 9.1, page 279):
           1. Compute forces at current configuration
           2. Compute acceleration from equation of motion
           3. Update velocity using central difference
           4. Update displacement
           """
           
           # Get state vectors
           u = globdat.state
           v = globdat.velo
           a = globdat.acce
           
           # Assemble forces at current configuration
           fint = assembleInternalForce(props, globdat)
           fext = assembleExternalForce(props, globdat)
           
           # Compute acceleration: a = M^{-1} (f_ext - f_int)
           residual = fext - fint
           a = residual / Mlumped
           
           # Apply boundary conditions
           a = applyConstraints(a)
           
           # Update velocity (central difference)
           v = v + a * dt
           
           # Update displacement
           u = u + v * dt
           
           # Update time
           time = time + dt
           
           # Check end time
           if time >= endTime:
               return False
           return True
   
       def checkStability():
           """Check CFL condition for stability.
           
           Time step must satisfy:
           dt < dt_crit = 2/ω_max
           
           where ω_max is the maximum natural frequency.
           """
           omega_max = sqrt(max(diag(K) / Mlumped))
           dt_crit = 2.0 / omega_max
           
           if dt > dt_crit:
               warning("Time step exceeds critical value")

**Full Implementation:**

See the complete implementation in:
`pyfem/solvers/ExplicitSolver.py <https://github.com/jjcremmers/PyFEM/blob/master/pyfem/solvers/ExplicitSolver.py>`_

Integration with PyFEM Framework
---------------------------------

Global Data Container
~~~~~~~~~~~~~~~~~~~~~

The ``globdat`` object contains all analysis data:

.. code-block:: python

   # State vectors
   globdat.state       # Displacement vector
   globdat.Dstate      # Displacement increment
   globdat.velo        # Velocity (dynamics)
   globdat.acce        # Acceleration (dynamics)
   
   # System data
   globdat.nodes       # Node coordinates
   globdat.dofs        # DOF manager
   globdat.elements    # Element groups
   
   # Solution control
   globdat.lam         # Load factor
   globdat.solverStatus # Status object with cycle, time, etc.

Assembly Routines
~~~~~~~~~~~~~~~~~

Use standard assembly functions from ``pyfem.fem.Assembly``:

.. code-block:: python

   from pyfem.fem.Assembly import (
       assembleTangentStiffness,  # Returns (K, fint)
       assembleInternalForce,      # Returns fint
       assembleExternalForce,      # Returns fext
       assembleMassMatrix          # Returns M
   )

Constraint Application
~~~~~~~~~~~~~~~~~~~~~~

Apply boundary conditions through the DOF manager:

.. code-block:: python

   # Constrain system matrix and vector
   K_c, f_c = globdat.dofs.constrain(K, f)
   
   # Constrain vector only
   a_c = globdat.dofs.constrainDofs(a)

Output Invocation
~~~~~~~~~~~~~~~~~

Output modules are called automatically, but you can control timing:

.. code-block:: python

   def writeOutput(self, globdat):
       """Called automatically by framework when solver yields."""
       pass

Registration and Usage
----------------------

File Location
~~~~~~~~~~~~~

Place solver in:

- ``pyfem/solvers/MySolver.py``

Import in __init__.py
~~~~~~~~~~~~~~~~~~~~~

Add to ``pyfem/solvers/__init__.py``:

.. code-block:: python

   from .MySolver import MySolver
   
   __all__ = [
       'MySolver',
       # ... other solvers
   ]

Using in Input Files
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   solver = 
   {
     type = "MySolver";
     
     # Solver-specific parameters
     maxCycle = 50;
     tol      = 1.0e-4;
   };

Testing and Validation
----------------------

Analytical Tests
~~~~~~~~~~~~~~~~

Verify against closed-form solutions:

.. code-block:: python

   def test_linear_problem(self):
       """Test solver on linear elastic problem."""
       # Should converge in one iteration
       pass

Convergence Studies
~~~~~~~~~~~~~~~~~~~

Test convergence rates for nonlinear problems:

1. **Quadratic convergence** for Newton-Raphson
2. **Superlinear convergence** for quasi-Newton methods
3. **Energy conservation** for dynamic solvers

Benchmark Problems
~~~~~~~~~~~~~~~~~~

Compare against standard benchmarks:

- Cook's membrane (plane stress)
- Cantilever beam (large displacements)
- Snap-through cylinder (arc-length)

Best Practices
--------------

1. **Check convergence**: Always monitor residual norms
2. **Handle exceptions**: Catch and report non-convergence
3. **Adaptive schemes**: Implement step size control when appropriate
4. **Logging**: Use the logger for user feedback
5. **Efficiency**: Reuse matrices when possible (e.g., mass matrix)
6. **Robustness**: Implement line search or damping for difficult problems
7. **Documentation**: Cite equations from the book

Common Pitfalls
---------------

1. **Forgetting to update load factor** between steps
2. **Not applying constraints** correctly
3. **Poor convergence criteria** (too loose or too tight)
4. **Incorrect time step** for explicit methods (CFL condition)
5. **Not resetting state** after non-convergence

References
----------

The theoretical foundation for solution algorithms can be found in:

*"Non-Linear Finite Element Analysis of Solids and Structures"*
by R. de Borst, M.A. Crisfield, J.J.C. Remmers and C.V. Verhoosel
John Wiley & Sons, 2012, ISBN 978-0470666449

Key chapters:

- Chapter 2: Nonlinear Finite Element Analysis
- Chapter 2.5: Path-Following Techniques (Arc-Length)
- Chapter 9: Dynamics and Time-Dependent Problems

See Also
--------

- :doc:`elements_dev` - Implementing element formulations
- :doc:`materials_dev` - Implementing material models
- :doc:`io_dev` - Implementing I/O modules
- :doc:`../solvers/overview` - Available solver types
