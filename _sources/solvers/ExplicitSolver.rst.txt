==============
ExplicitSolver
==============

The ``ExplicitSolver`` advances dynamics via a central-difference explicit
time-integration scheme using a lumped mass matrix. It updates velocity and
displacement, computes accelerations from the residual, and reports kinetic
energy.

Overview
--------

Solver type: ``ExplicitSolver``

- Method: explicit central differences with half-step velocity updates
- Mass: uses lumped mass matrix from assembly
- Loads: scalar load factor ``lam(t)`` applied to ``fhat``
- Termination: run stops at ``maxCycle``

Parameters
----------

Mandatory Parameters
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``type``
     - Must be set to ``"ExplicitSolver"``

Optional Parameters
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 75
   :header-rows: 1
   :width: 100%

   * - Parameter
     - Description
   * - ``dtime``
     - Time step size (required for stable integration)
   * - ``maxCycle``
     - Maximum cycles to run (default: unlimited)
   * - ``lam``
     - String expression for load factor vs time, e.g., ``"t"`` or ``"1.0"``

Algorithm Notes
---------------

- Updates: ``velo += 0.5*dt*acce``; ``disp += dt*velo``; recompute ``f_int``;
  compute accelerations from ``M_lumped acce = lam*fhat - f_int``; finalize
  ``velo`` with half-step.
- Sets constraint factor from current ``lam``.
- Prints kinetic energy periodically.

Examples
--------

- ``examples/ch05/StressWave_20x20.pro``

See Also
--------
- :doc:`../solvers`
- :doc:`../io/MeshWriter`
